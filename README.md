![Reedsolomon Banner](https://raw.githubusercontent.com/Archange-py/reedsolomon/refs/heads/main/resources/pictures/miniature_2.png)

Image generated with DALL-E

# Reed-Solomon Codes Implementation

[![GitHub Stars](https://img.shields.io/github/stars/Archange-py/reedsolomon.svg)](https://github.com/Archange-py/reedsolomon/stargazers)
[![GitHub Forks](https://img.shields.io/github/forks/Archange-py/reedsolomon.svg)](https://github.com/Archange-py/reedsolomon/network)
[![GitHub Issues](https://img.shields.io/github/issues/Archange-py/reedsolomon.svg)](https://github.com/Archange-py/reedsolomon/issues)
[![GitHub Pull Requests](https://img.shields.io/github/issues-pr/Archange-py/reedsolomon.svg)](https://github.com/Archange-py/reedsolomon/pulls)
[![GitHub License](https://img.shields.io/github/license/Archange-py/reedsolomon.svg)](https://github.com/Archange-py/reedsolomon/blob/main/LICENSE)
[![Last Commit](https://img.shields.io/github/last-commit/Archange-py/reedsolomon.svg)](https://github.com/Archange-py/reedsolomon/commits/main)
![Tests](https://github.com/Archange-py/reedsolomon/actions/workflows/python-tests.yml/badge.svg)
![Visiteurs](https://visitor-badge.laobi.icu/badge?page_id=Archange-py.reedsolomon)

[![PyPI Version](https://img.shields.io/pypi/v/reedsolomon.svg)](https://pypi.org/project/reedsolomon/)
[![PyPI - Downloads](https://img.shields.io/pypi/dm/reedsolomon.svg)](https://pypi.org/project/reedsolomon/)
[![Python Versions](https://img.shields.io/pypi/pyversions/reedsolomon.svg)](https://pypi.org/project/reedsolomon/)

Welcome to the **Reed-Solomon codes** GitHub repository! This project is dedicated to exploring Reed-Solomon codes. These codes are crucial in ensuring data integrity in digital communication and storage systems.

## Table of Contents
***
1. [General Info](#general-info)
2. [Getting Started](#getting-started)  
    - [Installation](#installation)  
    - [Dependencies](#dependencies)
    - [Usage](#usage)
3. [Repository Structure](#repository-structure)
4. [Links](#links)
5. [FAQs](#faqs)

## General Info

Reed-Solomon codes are a type of error-correcting code that operates on blocks of data. They are capable of detecting and correcting multiple symbol errors, making them invaluable in applications such as:

- QR Codes
- Data transmission in satellite and deep-space communication
- Storage devices (CDs, DVDs, Blu-rays)
- RAID storage systems

At present, only algorithms for polynomials and Galois fields are being developed. The construction of a mother class for the Reed-Solomon code will follow naturally.


## Getting Started

### Installation
***

You can use this command from PyPi:
- `pip install reedsolomon`


### Dependencies
***

Use this command in cmd if you are in a virtual environment:
- `pip install -r requirements.txt`


### Usage  
***

Here is examples to use the Polynomial file:  
- [`example_polynomial.ipynb`](./src/reedsolomon/examples/notebooks/examples_polynomial.ipynb)

The Galois Field file:
- [`example_galois_field.ipynb`](./src/reedsolomon/examples/notebooks/examples_galois_field.ipynb)

And the Reedsolomon file:
- [`example_reed_solomon.ipynb`](./src/reedsolomon/examples/notebooks/examples_reed_solomon.ipynb)


## Repository Structure

File composition in the repository :

```plaintext
.  
├── resources  
│   ├── documentation # all .pdf and .html  
│   └── pictures # for the README  
│
├── tests  
│   ├── test_reed_solomon.ipynb 
|   ├── test_galois_field.ipynb 
│   └── test_polynomial.ipynb  
|
├── src  
│   └── reedsolomon
│       ├── examples/notebooks  
│       │   ├── example_reed_solomon.ipynb 
|       │   ├── example_galois_field.ipynb 
│       │   └── example_polynomial.ipynb  
│       │
│       ├── __init__.py
│       ├── polynomial.py
│       ├── reedsolomon.py
│       └── galoisfield.py
│
├── MANIFEST.in
├── LICENSE.md
├── README.md
├── TODO.md
├── requirements.txt
├── pyproject.toml
└── .gitignore
```

## Links

Here you'll find the various resources used to complete this project:

- Resources sites
    - [Reed-Solomon Codes (Wikipedia)](https://en.wikipedia.org/wiki/Reed%E2%80%93Solomon_error_correction)

    - [Reed-Solomon Codes for coders (Wikiversity)](https://en.wikiversity.org/wiki/Reed%E2%80%93Solomon_codes_for_coders)

    - [An introduction to Reed-Solomon codes: principles, architecture and implementation](https://www.cs.cmu.edu/~guyb/realworld/reedsolomon/reed_solomon_codes.html)

    - [Finite Field Arithmetic and Reed Solomon Coding - Russ Cox](https://research.swtch.com/field)

- GitHub resources:
    - [Pure Python implementation](https://github.com/lrq3000/unireedsolomon)

    - [Optimized Python implementation](https://github.com/tomerfiliba-org/reedsolomon)

- [PDF resources](./Resource/Documentation/):
    - [CRC and Reed Solomon ECC - Jeff Reid.pdf](./Resource/Documentation/CRC_and_Reed_Solomon_ECC%20-%20Jeff%20Reid.pdf)

    - [Implementing Reed Solomon - Andrew Brown.pdf](./Resource/Documentation/Implementing_Reed_Solomon%20-%20Andrew%20Brown.pdf)

    - [Decoding Reed Solomon - Bruce Maggs.pdf](./Resource/Documentation/Decoding_Reed_Solomon%20-%20Bruce%20Maggs.pdf)

    - [Reed Solomon codes for coders - Wikiversity.pdf](./Resource/Documentation/Reed_Solomon_codes_for_coders%20-%20Wikiversity.pdf)

    - [Reed Solomon Scribe - Inconnu(e).pdf](./Resource/Documentation/Reed_Solomon_Scribe%20-%20Inconnu(e).pdf)

    - [WHP 031 - C.K.P. Clarke.pdf](./Resource/Documentation/WHP%20031%20-%20C.%20K.%20P.%20Clarke.pdf)

## FAQs

For the moment there is none.