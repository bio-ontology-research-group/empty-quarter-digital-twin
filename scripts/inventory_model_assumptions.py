"""Index numerical literals for review; inclusion is not evidence of calibration.

Run from any directory: uv run python scripts/inventory_model_assumptions.py
This inventories implementation constants, including structural/index constants.
Scientific interpretation and evidence status are described in Supplementary Note 4.
"""
import ast
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def inventory():
    files = sorted((ROOT / 'src/eq_biodt').rglob('*.py'))
    records = []
    hashes = {}
    for path in files:
        source = path.read_text()
        lines = source.splitlines()
        name = str(path.relative_to(ROOT))
        hashes[name] = hashlib.sha256(path.read_bytes()).hexdigest()
        tree = ast.parse(source)
        contexts = {}
        def walk(node, scope='module'):
            if isinstance(node, (ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)):
                scope += '.' + node.name
            if isinstance(node, ast.Constant) and type(node.value) in (int, float, complex):
                # Complex literals would need an explicit serialized representation.
                key = (node.lineno, scope)
                contexts.setdefault(key, []).append(repr(node.value))
            for child in ast.iter_child_nodes(node):
                walk(child, scope)
        walk(tree)
        for (line, scope), values in sorted(contexts.items()):
            records.append(dict(file=name, line=line, scope=scope,
                                source_line=lines[line-1].strip(), numeric_literals=values))
    return dict(description='All numerical literal occurrences in Python source. Includes indexing, dimensions, defaults, thresholds and coefficients. This is an implementation inventory, not a list of empirically estimated parameters. No evidence status is inferred from proximity to a citation.',
                coverage='src/eq_biodt/**/*.py; excludes data values, ontology axioms and browser code, which require their separate audits.',
                source_sha256=hashes, records=records)

if __name__ == '__main__':
    result = inventory()
    output = ROOT / 'reproducibility/outputs/model_coefficient_inventory.json'
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(f'{len(result["records"])} source locations in {len(result["source_sha256"])} files: {output}')
