# Changelog fragments

Every change lands here as one small file, and `towncrier` folds them into
`CHANGELOG.md` when a version is cut. The filename is `<issue-or-pr>.<type>.md`
(prefix the name with `+` when there is no issue, e.g. `+deprecation.added.md`), and the content is one sentence in the
past tense, the way it should read in the changelog.

| Type | Use it for |
| ---- | ---------- |
| `added` | a new public name or a new capability |
| `changed` | a change in behaviour that stays compatible |
| `deprecated` | something that now warns and will be removed |
| `removed` | something gone, only in a major version |
| `fixed` | a bug fixed |
| `security` | a security fix |
| `docs` | documentation only |
| `tests` | tests only |

```console
$ towncrier build --draft --version 1.0.0   # preview the section
$ towncrier build --version 1.0.0           # write it into CHANGELOG.md
```

`CHANGELOG.md` keeps its hand-written history above the version towncrier
generates; the fragments are deleted once they are folded in.
