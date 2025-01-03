# Building

## Conventional commits

This repo uses [conventional commits](https://www.conventionalcommits.org) formatted messages.

Basically, it means that your commit messages should be [atomic](https://en.wikipedia.org/wiki/Atomic_commit#Atomic_commit_convention)
and prepended with the type of change it introduces.

This is used for creating [CHANGELOG.md](CHANGELOG.md) and the changelog
object used by the JavaScript worker version.

For example

- `fix: Off by one error in output`
- `feat: Support multilayer files`
- `docs: Update installation instructions for Windows`

Other commit types are e.g. `perf`, `refactor`, `test` and `build`.

### Scopes

Scopes in conventional commits is appended to the commit type.

We use the scopes `js` and `python`, for example

- `fix(js): Check that error handler is a function`
- `docs(python): Update Python documentation`


## Releasing new versions

0. Commit everything with conventional commit messages
1. Install `standard-version` by running `npm i`
2. Run `npm run release -- --dry-run` to run [standard-version](https://github.com/conventional-changelog/standard-version)
    - Double check the output and that the `bumpFiles` has been updated
3. Run `npm run release` (Note: do not amend the release commit, the tag will point to the wrong commit!)
4. Run `git push --follow-tags`


### JavaScript

Building requires [Emscripten](https://emscripten.org/docs/getting_started/downloads.html).

In short:

```
git clone https://github.com/emscripten-core/emsdk.git
cd emsdk
./emsdk install latest
./emsdk activate latest
source ./emsdk_env.sh
```

We also need Infomap to extract the command line arguments for Infomap Online.

To build:

0. Follow the [release workflow](#releasing-new-versions) before building
1. Install deps with `npm i`
2. Run `make js-worker`
    - This creates `build/js/infomap.worker.js`
    - Copies the worker to `interfaces/js/src/worker`
    - Runs `npm run build` which bundles the worker with the js source files and copies them to `dist`
    - Copies the js README.md to the root
3. To test, run `make js-test`
    - Runs `npm pack` and extracts the `tgz` file.
    - Copies `packages/dist/index.js` to `examples/js/`.
    - Replaces the script source from the CDN to the local `./index.js`.
4. Run `npm publish`
    - Optionally run `make js-clean`


### Python

Building requires [Swig](http://swig.org) and [Sphinx](https://www.sphinx-doc.org)

On macOS, install with `brew install swig` and `brew install sphinx-doc`.

To build:

0. Follow the [release workflow](#releasing-new-versions) before building
1. Run `make python`
    - Runs `swig`
    - Creates `package-meta.py`
    - Copies python package files to `build/py`
2. To test the build, run `make py-test`
3. Test publish with `make pypitest-publish`
    - Install in a clean environment `pip3 --no-cache-dir install --index-url https://test.pypi.org/simple/ infomap`
4. Run `make pypi-dist`
5. Publish with `make pypi-publish`

Generate documentation:

0. Follow the [release workflow](#releasing-new-versions) before generating documentation
1. Run `make py-doc` which generates the documentation for Github pages. The front page is generated from `README.rst`
2. Commit the documentation with a `docs(python)` scoped commit.
