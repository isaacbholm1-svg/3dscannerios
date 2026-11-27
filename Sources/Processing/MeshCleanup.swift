import Foundation
import simd

/// Mesh clean-up and texture-label smoothing for LiDAR scans.
enum MeshCleanup {
    /// Merges vertices closer than `tolerance` (ARKit splits the mesh into
    /// anchors that each duplicate their border vertices) and drops triangles
    /// that collapse.
    static func weld(_ input: MeshData, tolerance: Float) -> MeshData {
        var out = MeshData()
        out.positions.reserveCapacity(input.positions.count)
        out.normals.reserveCapacity(input.positions.count)
        var remap = [UInt32](repeating: 0, count: input.positions.count)
        var cells: [SIMD3<Int32>: UInt32] = [:]
        cells.reserveCapacity(input.positions.count)
        let inv = 1 / tolerance
        for (i, p) in input.positions.enumerated() {
            let key = SIMD3<Int32>(Int32((p.x * inv).rounded()), Int32((p.y * inv).rounded()), Int32((p.z * inv).rounded()))
            if let existing = cells[key] {
                remap[i] = existing
                if i < input.normals.count { out.normals[Int(existing)] += input.normals[i] }
            } else {
                let index = UInt32(out.positions.count)
                cells[key] = index
                remap[i] = index
                out.positions.append(p)
                out.normals.append(i < input.normals.count ? input.normals[i] : SIMD3(0, 1, 0))
            }
        }
        for i in 0..<out.normals.count {
            let n = out.normals[i]
            let len = simd_length(n)
            out.normals[i] = len > 0 ? n / len : SIMD3(0, 1, 0)
        }
        let T = input.indices.count / 3
        out.indices.reserveCapacity(input.indices.count)
        out.faceClasses.reserveCapacity(T)
        for t in 0..<T {
            let a = remap[Int(input.indices[3 * t])]
            let b = remap[Int(input.indices[3 * t + 1])]
            let c = remap[Int(input.indices[3 * t + 2])]
            if a == b || b == c || a == c { continue }
            out.indices.append(a)
            out.indices.append(b)
            out.indices.append(c)
            out.faceClasses.append(t < input.faceClasses.count ? input.faceClasses[t] : 0)
        }
        return out
    }

    /// Compressed vertex -> neighbouring vertices lists.
    private static func vertexNeighbours(indices: [UInt32], vertexCount n: Int) -> (start: [Int], list: [UInt32]) {
        var count = [Int](repeating: 0, count: n + 1)
        let T = indices.count / 3
        for t in 0..<T {
            for j in 0..<3 { count[Int(indices[3 * t + j])] += 2 }
        }
        var start = [Int](repeating: 0, count: n + 1)
        for i in 0..<n { start[i + 1] = start[i] + count[i] }
        var fill = start
        var list = [UInt32](repeating: 0, count: start[n])
        for t in 0..<T {
            for j in 0..<3 {
                let v = Int(indices[3 * t + j])
                list[fill[v]] = indices[3 * t + (j + 1) % 3]
                list[fill[v] + 1] = indices[3 * t + (j + 2) % 3]
                fill[v] += 2
            }
        }
        return (start, list)
    }

    /// Taubin smoothing: removes LiDAR bumpiness on flat surfaces without
    /// shrinking the mesh, then recomputes vertex normals.
    static func taubinSmooth(_ mesh: inout MeshData, iterations: Int, lambda: Float = 0.5, mu: Float = -0.53) {
        let n = mesh.positions.count
        guard n > 0, mesh.indices.count >= 3 else { return }
        let (start, list) = vertexNeighbours(indices: mesh.indices, vertexCount: n)
        var next = mesh.positions
        for _ in 0..<iterations {
            for factor in [lambda, mu] {
                for v in 0..<n {
                    let s = start[v], e = start[v + 1]
                    guard e > s else { next[v] = mesh.positions[v]; continue }
                    var sum = SIMD3<Float>(repeating: 0)
                    for i in s..<e { sum += mesh.positions[Int(list[i])] }
                    let avg = sum / Float(e - s)
                    next[v] = mesh.positions[v] + factor * (avg - mesh.positions[v])
                }
                swap(&mesh.positions, &next)
            }
        }
        // Fresh area-weighted vertex normals.
        var normals = [SIMD3<Float>](repeating: .zero, count: n)
        let T = mesh.indices.count / 3
        for t in 0..<T {
            let a = Int(mesh.indices[3 * t]), b = Int(mesh.indices[3 * t + 1]), c = Int(mesh.indices[3 * t + 2])
            var fn = simd_cross(mesh.positions[b] - mesh.positions[a], mesh.positions[c] - mesh.positions[a])
            // Keep ARKit's outward orientation.
            if simd_dot(fn, mesh.normals[a] + mesh.normals[b] + mesh.normals[c]) < 0 { fn = -fn }
            normals[a] += fn
            normals[b] += fn
            normals[c] += fn
        }
        for v in 0..<n {
            let len = simd_length(normals[v])
            mesh.normals[v] = len > 0 ? normals[v] / len : mesh.normals[v]
        }
    }

    /// Chooses one photo per triangle from its top-k candidates, trading a bit
    /// of per-triangle quality for agreement with neighbouring triangles.
    static func smoothLabels(indices: [UInt32], vertexCount n: Int, candFrame: [Int32], candScore: [Float],
                             k: Int, iterations: Int = 4, neighbourWeight: Float = 0.8) -> [Int32] {
        let T = indices.count / 3
        var label = [Int32](repeating: -1, count: T)
        for t in 0..<T { label[t] = candFrame[t * k] }
        guard T > 1 else { return label }

        // Vertex -> triangles, then triangle neighbours through shared vertices.
        var vcount = [Int](repeating: 0, count: n + 1)
        for i in indices { vcount[Int(i)] += 1 }
        var vstart = [Int](repeating: 0, count: n + 1)
        for v in 0..<n { vstart[v + 1] = vstart[v] + vcount[v] }
        var vfill = vstart
        var vtris = [Int32](repeating: 0, count: vstart[n])
        for t in 0..<T {
            for j in 0..<3 {
                let v = Int(indices[3 * t + j])
                vtris[vfill[v]] = Int32(t)
                vfill[v] += 1
            }
        }

        for _ in 0..<iterations {
            var changed = 0
            for t in 0..<T {
                let base = t * k
                let best = candScore[base]
                guard best > 0, candFrame[base + 1] >= 0 else { continue }
                // How many neighbouring triangles already use each candidate photo.
                var votes = [Int](repeating: 0, count: k)
                var total = 0
                for j in 0..<3 {
                    let v = Int(indices[3 * t + j])
                    for i in vstart[v]..<vstart[v + 1] {
                        let u = Int(vtris[i])
                        if u == t { continue }
                        total += 1
                        let l = label[u]
                        for c in 0..<k where candFrame[base + c] == l { votes[c] += 1 }
                    }
                }
                guard total > 0 else { continue }
                var choice = label[t]
                var choiceValue: Float = -1
                for c in 0..<k {
                    let f = candFrame[base + c]
                    if f < 0 { break }
                    let quality = candScore[base + c] / best
                    if quality < 0.35 { continue }
                    let value = quality + neighbourWeight * Float(votes[c]) / Float(total)
                    if value > choiceValue {
                        choiceValue = value
                        choice = f
                    }
                }
                if choice != label[t] {
                    label[t] = choice
                    changed += 1
                }
            }
            if changed == 0 { break }
        }
        return label
    }
}
