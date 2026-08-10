# NMR Restraints Mapping

This document maps the experimental NMR restraints to the corresponding atom names in the provided structure files (`cx4.mol2` and `atto655.mol2`).

## 1. sCX4 (Host) Ha Atoms

The restraints are measured between the Ha atoms of the sCX4 host and specific H atoms of the Atto655 dye.
There are 8 Ha atoms across the 4 subunits of sCX4.

| NMR Naming | Subunit | cx4.mol2 Atom Name |
| :--- | :---: | :--- |
| Ha | 1 | `/cx4///CX4`1/H61` |
| Ha | 1 | `/cx4///CX4`1/H57` |
| Ha | 2 | `/cx4///CX4`1/H58` |
| Ha | 2 | `/cx4///CX4`1/H62` |
| Ha | 3 | `/cx4///CX4`1/H60` |
| Ha | 3 | `/cx4///CX4`1/H64` |
| Ha | 4 | `/cx4///CX4`1/H63` |
| Ha | 4 | `/cx4///CX4`1/H59` |

<br>

## 2. Atto655 (Guest) H Atoms

The following table maps the Atto655 NMR naming scheme to the specific heavy atoms (C) and their attached protons (H1, H2, H3).

| NMR Naming | Atto655.mol2 (Heavy Atom) | H1 (mol2) | H2 (mol2) | H3 (mol2) |
| :---: | :--- | :--- | :--- | :--- |
| 1 | `/atto655///DDS`1/C19` | | | |
| 2 | `/atto655///DDS`1/C18` | `/atto655///DDS`1/H17` | `/atto655///DDS`1/H16` | |
| 3 | `/atto655///DDS`1/C17` | `/atto655///DDS`1/H14` | `/atto655///DDS`1/H15` | |
| 4 | `/atto655///DDS`1/C16` | `/atto655///DDS`1/H12` | `/atto655///DDS`1/H13` | |
| 5 | `/atto655///DDS`1/C10` | | | |
| 6 | `/atto655///DDS`1/C9` | `/atto655///DDS`1/H7` | | |
| 7 | `/atto655///DDS`1/C8` | | | |
| 8 | `/atto655///DDS`1/C5` | | | |
| 9 | `/atto655///DDS`1/C4` | `/atto655///DDS`1/H6` | | |
| 10 | `/atto655///DDS`1/C3` | | | |
| 11 | `/atto655///DDS`1/C26` | | | |
| 12 | `/atto655///DDS`1/C25` | `/atto655///DDS`1/H33` | | |
| 13 | `/atto655///DDS`1/C6` | | | |
| 14 | `/atto655///DDS`1/C7` | | | |
| 15 | `/atto655///DDS`1/C12` | `/atto655///DDS`1/H8` | | |
| 16 | `/atto655///DDS`1/C11` | | | |
| 17 | `/atto655///DDS`1/C13` | `/atto655///DDS`1/H9` | | |
| 18 | `/atto655///DDS`1/C14` | `/atto655///DDS`1/H10` | `/atto655///DDS`1/H11` | |
| 19 | `/atto655///DDS`1/C15` | | | |
| 20 | `/atto655///DDS`1/C22` | `/atto655///DDS`1/H25` | `/atto655///DDS`1/H27` | `/atto655///DDS`1/H26` |
| 21 | `/atto655///DDS`1/C23` | `/atto655///DDS`1/H30` | `/atto655///DDS`1/H28` | `/atto655///DDS`1/H29` |
| 22 | `/atto655///DDS`1/C24` | `/atto655///DDS`1/H31` | `/atto655///DDS`1/H32` | |
| 23 | `/atto655///DDS`1/C27` | `/atto655///DDS`1/H34` | `/atto655///DDS`1/H35` | |
| 24 | `/atto655///DDS`1/C28` | `/atto655///DDS`1/H37` | `/atto655///DDS`1/H36` | |
| 25 | `/atto655///DDS`1/C29` | `/atto655///DDS`1/H39` | `/atto655///DDS`1/H38` | |
| 26 | `/atto655///DDS`1/C2` | `/atto655///DDS`1/H5` | `/atto655///DDS`1/H4` | |
| 27 | `/atto655///DDS`1/C1` | `/atto655///DDS`1/H3` | `/atto655///DDS`1/H1` | `/atto655///DDS`1/H2` |

<br>

## 3. NOESY Restraints (Attractive)

These are H atoms on the Atto655 dye that are close to the Ha atoms of sCX4. The defined upper distance bound is specified.

| Atto655 NMR Naming | Atto655 Protons (mol2) | Target sCX4 Atoms | Peak Class | Upper Distance (Å) |
| :---: | :--- | :--- | :--- | :---: |
| **12** | `/atto655///DDS`1/H33` | 8 $\times$ Ha atoms | Strong | 3.5 |
| **6** | `/atto655///DDS`1/H7` | 8 $\times$ Ha atoms | Weak | 3.55 |
| **9** | `/atto655///DDS`1/H6` | 8 $\times$ Ha atoms | Medium | 3.5 |
| **27** | `/atto655///DDS`1/H3`, `/atto655///DDS`1/H1`, `/atto655///DDS`1/H2` | 8 $\times$ Ha atoms | Strong | 3.5 |
| **24** | `/atto655///DDS`1/H37`, `/atto655///DDS`1/H36` | 8 $\times$ Ha atoms | Weak | 5.0 |
| **23** | `/atto655///DDS`1/H34`, `/atto655///DDS`1/H35` | 8 $\times$ Ha atoms | Medium | 5.0 |
| **25** | `/atto655///DDS`1/H39`, `/atto655///DDS`1/H38` | 8 $\times$ Ha atoms | Medium | 5.0 |
| **26** | `/atto655///DDS`1/H5`, `/atto655///DDS`1/H4` | 8 $\times$ Ha atoms | Medium | 5.0 |

<br>

## 4. CSP Restraints (Repulsive)

These are specific H atoms of the Atto655 dye that are *not* in contact with the Ha atoms of the sCX4 cage, and thus act as repulsive restraints.

| Atto655 NMR Naming | Atto655 Protons (mol2) | Rule |
| :---: | :--- | :--- |
| **3** | `/atto655///DDS`1/H14`, `/atto655///DDS`1/H15` | Repulsion from Ha atoms |
| **4** (4-1, 4-2) | `/atto655///DDS`1/H12`, `/atto655///DDS`1/H13` | Repulsion from Ha atoms |
| **6** | `/atto655///DDS`1/H7` | Repulsion from Ha atoms |
| **18** (18-1, 18-2) | `/atto655///DDS`1/H10`, `/atto655///DDS`1/H11` | Repulsion from Ha atoms |
| **20** | `/atto655///DDS`1/H25`, `/atto655///DDS`1/H27`, `/atto655///DDS`1/H26` | Repulsion from Ha atoms |
| **21** | `/atto655///DDS`1/H30`, `/atto655///DDS`1/H28`, `/atto655///DDS`1/H29` | Repulsion from Ha atoms |
