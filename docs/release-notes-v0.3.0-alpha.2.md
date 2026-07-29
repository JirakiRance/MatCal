# MatCal v0.3.0-alpha.2 Release Notes

`v0.3.0-alpha.2` is the next release candidate after `v0.3.0-alpha.1`.

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
- Records the missing repository license as a release blocker in `docs/licensing-decision.md`.

## CI Status Requirement

Do not create `v0.3.0-alpha.2` until a pushed alpha.2 candidate commit has green GitHub Linux GCC and Windows MSVC CI.

If Linux CI still fails, the workflow should provide:

- verbose rerun output for failed tests;
- GitHub annotations containing the failed CTest names and assertion summaries;
- `linux-debug-ctest-logs` artifact for Debug failures;
- `linux-release-ctest-logs` artifact for Release failures when Release tests run.

## SFL Upgrade Note

SFL is temporarily using `v0.3.0-alpha.1`. Once `v0.3.0-alpha.2` is green, SFL should update only the MatCal gitlink/tag to the new candidate and continue linking the required MatCal target, primarily `MatCal::Linalg` for Skyline/LDLT use. SFL must not copy MatCal sources or add MatCal as an internal subtree.

After updating the gitlink, SFL should rerun native TRI3/QUAD4, Result, and course model regressions.

## Known Release Blockers

- Linux GCC CI failure root cause is not yet known because full GitHub job logs were not accessible through the unauthenticated/admin-limited API.
- Repository license is not declared. See `docs/licensing-decision.md`.
