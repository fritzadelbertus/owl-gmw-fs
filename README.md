# OWL-GMW-FS

Status: Research Prototype

Implementation of the OWL signature scheme using the
Goldreich–Micali–Wigderson identification protocol transformed
via the Fiat–Shamir heuristic.

## Overview

OWL-GMW-FS is an experimental post-quantum digital signature scheme constructed from a group-action-based identification protocol.

The repository contains

- key generation
- signing
- verification
- unit tests
- benchmarking
- customization

This implementation accompanies our research on constructing group-action-based signatures inspired by ALTEQ and HAWK.

## Security Notice

This project is an experimental research implementation.

It has not been independently audited and should not be used in production environments.

## Features

- General GMW-FS Framework
- OWL-mLIP: GMW-FS instantiated using mLIGA (group action derived from module Lattice Isomorphism Problem)
- OWL-LIP: GMW-FS instantiated using LIGA (group action derived from Lattice Isomorphism Problem)
- Unit testing to test any GMW-FS based digital signature protocols
- Benchmarking

## Installation

Below are the steps run in terminal to install and setup this project

```bash
git clone https://github.com/yourname/OWL-GMW-FS.git
cd OWL-GMW-FS

python -m venv .venv
source .venv/bin/activate

pip install -r requirements.txt
```

## How to use: OWL-mLIP

### Generating Key, Signing Messages, Verifying Signatures

This code below shows how to use the three digital signature algorithms of OWL-mLIP

```python
from owl_mLIP.main import owl_mLIP 

pk, sk, orb = owl_mLIP.keygen()

msg = os.urandom(32) # byte object

sig = owl_mLIP.sign(sk, orb, msg)

result = owl_mLIP.verify(pk, orb, msg, sig)
```

### Customizing Parameters

OWL-mLIP have several modifiable parameters in ```owl-mLIP/params.py```

- ```LAMBDA```: n-bit security parameter
- ```C```: number of group & set elements sampled during key generation 
- ```K```: parameter for unbalanced challenge
- ```R```: number of group elements sample during signing
- ```LOGN```: parameter derived from HAWK protocol

## How to use: OWL-LIP

### Generating Key, Signing Messages, Verifying Signatures

This code below shows how to use the three digital signature algorithms of OWL-LIP

```python
from owl_LIP.main import owl_LIP
from owl_LIP.params import orbit

pk, sk, orb = owl_LIP.keygen(orbit)

msg = os.urandom(32) # byte object

sig = owl_LIP.sign(sk, orb, msg)

result = owl_LIP.verify(pk, orb, msg, sig)
```

### Customizing Parameters

OWL-LIP have several modifiable parameters in ```owl-LIP/params.py```

- ```LAMBDA```: n-bit security parameter
- ```C```: number of group & set elements sampled during key generation 
- ```K```: parameter for unbalanced challenge
- ```R```: number of group elements sample during signing
- ```L```: number of walks to generate the random matrix in ```sample.py```

## How to create you own GMW-FS based digital signature protocol

If you have a cryptographic group action, a custom gmwfs can easily be instantiated. Create a copy of ```custom_gmwfs``` directory and fill in every each files.

## Benchmarks

The benchmark result for this project used a machine describe as follows 

CPU: AMD Ryzen 7 7840HS (8 cores, 16 threads)
RAM: 32 GB DDR5
Operating System: Arch Linux
Python: 3.14
NumPy: 2.3.1
BLAS: OpenBLAS


To run basic benchmarks:
```bash
python -m benchmarks.benchmark --scheme (owl_lip/owl_mlip)
```

Benchmarks are currently under development and will be included in a future release.

## Testing

Files inside ```tests``` used to verify the correctness of each cryptographic group action and the gmw-fs components

To run test do
```bash
python -m pytest -v tests/test_name.py
```

## Citation
HAWK and ALTEQ in Group Action Based GMW-FS Design

(BibTeX will be added upon publication)

## Reference
- Hawk: Module LIP make Lattice Signature Fast, Compact, and Simple
- ALTEQ Specification

## License
TBD


