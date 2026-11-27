# 3D Scanner for iOS

A SwiftUI app for capturing 3D scans on iPhone and iPad: objects, LiDAR spaces, rooms (RoomPlan) and photo sets, with a scan library, a model viewer and exporters.

## Building

There is no checked-in Xcode project. [XcodeGen](https://github.com/yonaskolb/XcodeGen) generates it from `project.yml`:

```sh
brew install xcodegen
xcodegen generate
open Scanner3D.xcodeproj
```

GitHub Actions (`.github/workflows/ios.yml`) builds an unsigned `Scanner3D.ipa` on every push. Download it from the run's Artifacts (or from the `scanner-latest` release for builds on `main`), then sign and sideload it with KSign, AltStore or Sideloadly.

Requires iOS 17+. Object, Space and Room modes need a LiDAR device (iPhone/iPad Pro).
