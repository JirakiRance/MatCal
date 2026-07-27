# Compatibility Policy

MatCal M0 establishes a compatibility baseline. It does not declare production readiness or ABI stability.

## Compatibility Levels

Behavior compatibility means existing documented legacy behavior keeps working unless it is a confirmed bug with an explicit migration path.

Source/API compatibility means existing public headers, namespaces, class names, function names, constructor forms, and argument order remain available in the current major line.

Binary/ABI compatibility means already compiled consumers can link against a new binary without recompilation. MatCal does not promise this today. The current library is effectively header-only, compiler-specific, and not ABI-governed.

## M0 Commitments

- Do not delete public headers.
- Do not rename public namespaces or classes.
- Do not remove old functions directly.
- Prefer adding safe new APIs and keeping deprecated wrappers for public mistakes.
- Keep `MatCal::Legacy` as the target for the current public surface.
- Do not overwrite or replace historical compiled binaries used by PT or other projects.
- Do not promise cross-compiler or cross-version ABI stability.

## Deprecation Strategy

When an existing public API is unsafe but likely used:

1. Add a replacement with clearer ownership, constness, status reporting, and numeric contract.
2. Keep the legacy entry point forwarding to the replacement where feasible.
3. Mark the legacy API deprecated only after tests and migration notes exist.
4. Remove only at a declared major-version boundary.

## Version Boundary

Suggested boundaries:

- `0.x`: compatibility baseline, tests, internal refactor, no ABI promise.
- `1.x`: formal source-compatible stable legacy facade plus selected new APIs.
- `2.0`: first possible removal window for deprecated APIs, with migration guide.
