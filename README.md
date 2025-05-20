# Path-based fingerprint

Path-based fingerprints calculated with the RDKit package Chem.RDKFingerprint. It is inspired in the Daylight fingerprint. As explained in the RDKit Book, the fingerprinting algorithm identifies all subgraphs in the molecule within a particular range of sizes, hashes each subgraph to generate a raw bit ID, mods that raw bit ID to fit in the assigned fingerprint size, and then sets the corresponding bit. 

This model was incorporated on 2021-09-17.

## Information
### Identifiers
- **Ersilia Identifier:** `eos7jio`
- **Slug:** `rdkit-fingerprint`

### Domain
- **Task:** `Representation`
- **Subtask:** `Featurization`
- **Biomedical Area:** `Any`
- **Target Organism:** `Not Applicable`
- **Tags:** `Fingerprint`, `Descriptor`

### Input
- **Input:** `Compound`
- **Input Dimension:** `1`

### Output
- **Output Dimension:** `2048`
- **Output Consistency:** `Fixed`
- **Interpretation:** Vector representation of small molecules

Below are the **Output Columns** of the model:
| Name | Type | Direction | Description |
|------|------|-----------|-------------|
| dimension_0000 | integer |  | RDKit fingerprint dimension 0 with 2048 bits and path length of 1-7 |
| dimension_0001 | integer |  | RDKit fingerprint dimension 1 with 2048 bits and path length of 1-7 |
| dimension_0002 | integer |  | RDKit fingerprint dimension 2 with 2048 bits and path length of 1-7 |
| dimension_0003 | integer |  | RDKit fingerprint dimension 3 with 2048 bits and path length of 1-7 |
| dimension_0004 | integer |  | RDKit fingerprint dimension 4 with 2048 bits and path length of 1-7 |
| dimension_0005 | integer |  | RDKit fingerprint dimension 5 with 2048 bits and path length of 1-7 |
| dimension_0006 | integer |  | RDKit fingerprint dimension 6 with 2048 bits and path length of 1-7 |
| dimension_0007 | integer |  | RDKit fingerprint dimension 7 with 2048 bits and path length of 1-7 |
| dimension_0008 | integer |  | RDKit fingerprint dimension 8 with 2048 bits and path length of 1-7 |
| dimension_0009 | integer |  | RDKit fingerprint dimension 9 with 2048 bits and path length of 1-7 |

_10 of 2048 columns are shown_
### Source and Deployment
- **Source:** `Local`
- **Source Type:** `External`
- **DockerHub**: [https://hub.docker.com/r/ersiliaos/eos7jio](https://hub.docker.com/r/ersiliaos/eos7jio)
- **Docker Architecture:** `AMD64`, `ARM64`
- **S3 Storage**: [https://ersilia-models-zipped.s3.eu-central-1.amazonaws.com/eos7jio.zip](https://ersilia-models-zipped.s3.eu-central-1.amazonaws.com/eos7jio.zip)

### Resource Consumption


### References
- **Source Code**: [https://github.com/rdkit/rdkit](https://github.com/rdkit/rdkit)
- **Publication**: [https://www.rdkit.org/docs/RDKit_Book.html](https://www.rdkit.org/docs/RDKit_Book.html)
- **Publication Type:** `Other`
- **Publication Year:** `2010`
- **Ersilia Contributor:** [miquelduranfrigola](https://github.com/miquelduranfrigola)

### License
This package is licensed under a [GPL-3.0](https://github.com/ersilia-os/ersilia/blob/master/LICENSE) license. The model contained within this package is licensed under a [BSD-3-Clause](LICENSE) license.

**Notice**: Ersilia grants access to models _as is_, directly from the original authors, please refer to the original code repository and/or publication if you use the model in your research.


## Use
To use this model locally, you need to have the [Ersilia CLI](https://github.com/ersilia-os/ersilia) installed.
The model can be **fetched** using the following command:
```bash
# fetch model from the Ersilia Model Hub
ersilia fetch eos7jio
```
Then, you can **serve**, **run** and **close** the model as follows:
```bash
# serve the model
ersilia serve eos7jio
# generate an example file
ersilia example -n 3 -f my_input.csv
# run the model
ersilia run -i my_input.csv -o my_output.csv
# close the model
ersilia close
```

## About Ersilia
The [Ersilia Open Source Initiative](https://ersilia.io) is a tech non-profit organization fueling sustainable research in the Global South.
Please [cite](https://github.com/ersilia-os/ersilia/blob/master/CITATION.cff) the Ersilia Model Hub if you've found this model to be useful. Always [let us know](https://github.com/ersilia-os/ersilia/issues) if you experience any issues while trying to run it.
If you want to contribute to our mission, consider [donating](https://www.ersilia.io/donate) to Ersilia!
