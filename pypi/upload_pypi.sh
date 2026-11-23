#!/bin/bash
python -m pip install build twine auditwheel patchelf

out=$(mktemp -d)
python -m build --outdir "$out/dist"
auditwheel repair "$out"/dist/*.whl -w "$out/wheelhouse"

python -m twine check "$out"/dist/*.tar.gz "$out"/wheelhouse/*.whl

python -m twine upload --repository pypi -u __token__ \
"$out"/dist/*.tar.gz \
"$out"/wheelhouse/*.whl
