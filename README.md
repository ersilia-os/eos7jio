# Path-based fingerprint

Path-based fingerprints calculated with the RDKit package Chem.RDKFingerprint. It is inspired in the Daylight fingerprint. As explained in the RDKit Book, the fingerprinting algorithm identifies all subgraphs in the molecule within a particular range of sizes, hashes each subgraph to generate a raw bit ID, mods that raw bit ID to fit in the assigned fingerprint size, and then sets the corresponding bit. 

## Identifiers

* EOS model ID: `eos7jio`
* Slug: `rdkit-fingerprint`

## Characteristics

* Input: `Compound`
* Input Shape: `Single`
* Task: `Representation`
* Output: `Descriptor`
* Output Type: `Float`
* Output Shape: `List`
* Interpretation: Vector representation of small molecules

## References

* [Publication](https://www.rdkit.org/docs/RDKit_Book.html)
* [Source Code](https://github.com/rdkit/rdkit)
* Ersilia contributor: [miquelduranfrigola](https://github.com/miquelduranfrigola)

## Citation

If you use this model, please cite the [original authors](https://www.rdkit.org/docs/RDKit_Book.html) of the model and the [Ersilia Model Hub](https://github.com/ersilia-os/ersilia/blob/master/CITATION.cff).

## License

This package is licensed under a GPL-3.0 license. The model contained within this package is licensed under a BSD-3.0 license.

Notice: Ersilia grants access to these models 'as is' provided by the original authors, please refer to the original code repository and/or publication if you use the model in your research.

## About Us

The [Ersilia Open Source Initiative](https://ersilia.io) is a Non Profit Organization ([1192266](https://register-of-charities.charitycommission.gov.uk/charity-search/-/charity-details/5170657/full-print)) with the mission is to equip labs, universities and clinics in LMIC with AI/ML tools for infectious disease research.

[Help us](https://www.ersilia.io/donate) achieve our mission!