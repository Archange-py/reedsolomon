<img src=".\Resource\Pictures\picture_2.png" alt="DALL-E Picture">  

Image generated with DALL-E

# Reed-Solomon Codes Implementation

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

There is none for the moment.


### Dependencies
***

Use this command in cmd:
- `pip install -r requirements.txt`

### Usage  
***

Here is examples to use the `polynomial.py` file:  
- [`example_polynomial.ipynb`](./Source/Examples/examples_polynomial.ipynb)

The galois field file:
- [`example_galois_field.ipynb`](./Source/Examples/examples_galois_field.ipynb)

And the reed solomon file:
- [`example_reed_solomon.ipynb`](./Source/Examples/examples_reed_solomon.ipynb)


## Repository Structure

File composition in the repository :

```plaintext
.  
├── Ressource  
│   ├── Documentation # all .pdf and .html  
│   └── Picture # for the README  
│
├── Source  
│   ├── Examples  
│   │   ├── example_reed_solomon.ipynb
│   │   ├── example_galois_field.ipynb
│   │   └── example_polynomial.ipynb
│   │
│   └── Tests
│       ├── test_reed_solomon.py
│       ├── test_galois_field.py
│       └── test_polynomial.py
│
├── README.md
├── LICENSE.md
├── TODO.md
├── notes.txt
├── requirements.txt
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