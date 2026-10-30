"""Install the checksum-pinned Linux x86_64 reasoner used for publication."""
from pathlib import Path
import hashlib
import platform
import urllib.request

# Keep in sync with KM_VERSION/KM_SHA256 in src/eq_biodt/v03/ontology.py (runtime digest check; tests assert equality).
VERSION = '1.5.3'
SHA256 = '2470b0415f1a78533024e15b5c317f9e743a98916b597739db93b9452a8dd8d4'
URL = f'https://github.com/bio-ontology-research-group/kobayashi-marust/releases/download/v{VERSION}/km-v{VERSION}-linux-x86_64'

def main():
    if platform.system() != 'Linux' or platform.machine() not in ('x86_64', 'AMD64'):
        raise SystemExit('Publication reproduction uses Linux x86_64. Use the Linux CI/container, or build km v1.5.3 and set KM_BIN; ordinary tool use needs no reasoner binary.')
    dest = Path(__file__).resolve().parents[1] / '.tools' / 'km'
    dest.parent.mkdir(exist_ok=True)
    data = dest.read_bytes() if dest.exists() else urllib.request.urlopen(URL, timeout=120).read()
    if hashlib.sha256(data).hexdigest() != SHA256:
        raise SystemExit('Reasoner checksum mismatch; refusing to install.')
    dest.write_bytes(data)
    dest.chmod(0o755)
    print(f'km {VERSION}: {dest} (SHA256 {SHA256})')

if __name__ == '__main__':
    main()
