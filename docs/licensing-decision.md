# Licensing Decision Required

MatCal currently has no confirmed repository license.

## Checks Performed

- Repository root was checked for `LICENSE`, `COPYING`, `NOTICE`, and `THIRD_PARTY` files; none are present.
- `README.md` and public source/header introductions were searched for license, copyright, MIT, Apache, BSD, GPL, and related terms; no license declaration was found.
- Git history was searched for license-like filenames and text. No project license declaration was found. Historical build artifact paths in old commits contain words such as `ContinuousSubmit`, but they are not license metadata.
- GitHub repository license endpoint was queried and returned `404 Not Found`, meaning GitHub does not detect a license file for the repository.
- No `.gitmodules` file is present.

## Decision Needed

The repository owner must choose and add a license before MatCal is treated as a redistributable public release package.

Codex must not infer or select a license on behalf of the owner. Common choices such as MIT, BSD, Apache-2.0, GPL, or proprietary/no-license all have different legal consequences.

## Third-Party Notice

No bundled third-party source code was identified in the current source tree during this pass. Do not create a `THIRD_PARTY` or `NOTICE` file that invents dependencies. If future vendored code or copied algorithms are added, record their licenses explicitly.

## Release Impact

`v0.3.0-alpha.2` can be prepared as a technical candidate, but publication or broader redistribution should remain blocked until the owner records the intended license.
