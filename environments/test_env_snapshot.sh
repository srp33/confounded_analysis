echo "=== PATH ==="
echo "$PATH" | tr ':' '\n' | head -n 10

echo "=== VIRTUAL_ENV ==="
echo "$VIRTUAL_ENV"

echo "=== PYTHONHOME ==="
echo "${PYTHONHOME:-<unset>}"

echo "=== command -v python ==="
command -v python

echo "=== command -v pip ==="
command -v pip

echo "=== command -v R ==="
command -v R

echo "=== which python, pip, R ==="
which python
which pip
which R

echo "=== env | grep -E 'VIRTUAL|PYTHON' ==="
env | grep -E 'VIRTUAL|PYTHON'
