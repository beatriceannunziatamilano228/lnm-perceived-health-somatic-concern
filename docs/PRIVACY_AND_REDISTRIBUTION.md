# Privacy and redistribution audit

The release intentionally excludes:

- participant identifiers;
- internal row-to-identifier matching files;
- participant-level clinical spreadsheets;
- individual lesion masks;
- participant-level lesion-connectivity matrices;
- normative-connectome subject data;
- raw SimNIBS logs and session files containing personal filesystem paths.

Included NIfTI files are group-level statistical or target maps. Included NumPy arrays are group maps or permutation null distributions, not participant records.
