# MatCal v0.3.0-alpha.2 Release Notes

`v0.3.0-alpha.2` is the released alpha.2 snapshot after `v0.3.0-alpha.1`.

`v0.3.0-alpha.1` points to `31f493789a5b2f191aef59c828ec7562a7f9412f` and has already been pushed. It must not be moved, overwritten, or deleted. Its GitHub Windows MSVC job passed, but its Linux GCC job failed in the Debug CTest step.

## API Compatibility

The public source API is intended to remain compatible with `v0.3.0-alpha.1`.

No new numerical algorithms, CMake targets, public classes, solver statuses, or solver options are added in this candidate. The CMake project version remains numeric `0.3.0`; the prerelease identity is expressed by the Git tag and release notes.

## Changes Since alpha.1

- Adds the alpha API manifest to the main branch baseline.
- Adds Linux CI diagnostics so failed CTest runs rerun failed tests verbosely.
- Emits failed Linux CTest names and `FAIL:` assertion summaries as GitHub annotations so they remain visible even when artifact downloads require repository permissions.
- Uploads Linux CTest `LastTest.log` and `LastTestsFailed.log` artifacts on failure when those files exist.
- Restricts Debug and Release rerun/upload diagnostics to the matching failed CTest step, avoiding spurious Release diagnostics after an earlier Debug failure.
- Fixes a Picard-Aitken cancellation edge case that made `matcal_roots_tests` fail on Linux GCC Debug: near the fixed point, the Aitken denominator is now treated as unstable when it is below the requested iterate tolerance, and the step falls back to the ordinary Picard update instead of amplifying roundoff.
- Records the then-missing repository license as a release blocker for the
  original alpha.2 candidate. The D1 maintenance cleanup later resolved the
  project license on `main` by adding the MIT `LICENSE` file.

## Published CI Status

The published tag points to `c1daae9d8785b82fd53370a7c903d6763a21dfbd`.
Both main and tag CI completed successfully:

- Main run: <https://github.com/JirakiRance/MatCal/actions/runs/30686518196>
- Tag run: <https://github.com/JirakiRance/MatCal/actions/runs/30686631859>
- Linux GCC: success
- Windows MSVC: success

## SFL Upgrade Note

SFL can update only the MatCal gitlink/tag to `v0.3.0-alpha.2` and continue
linking the required MatCal target, primarily `MatCal::Linalg` for Skyline/LDLT
use. SFL must not copy MatCal sources or add MatCal as an internal subtree.

After updating the gitlink, SFL should rerun native TRI3/QUAD4, Result, and course model regressions.

## Release Asset Notes

- The published `v0.3.0-alpha.2` tag points to
  `c1daae9d8785b82fd53370a7c903d6763a21dfbd`, which predates the D1
  maintenance cleanup. That tag source snapshot does not contain the later MIT
  `LICENSE` file.
- Current `main` contains the MIT license. Do not move the published tag to add
  post-release files.
