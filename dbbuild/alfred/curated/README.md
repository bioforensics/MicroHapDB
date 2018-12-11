The microhaplotype score data in this directory are "curated" in the sense that it was never provided in a convenient machine-readable format by the authors.
It had to be copy-and-pasted from the article HTML and manually re-formatted, which is always an error-prone process.

The microhap locus coordinates here were derived from the GRCh38 coordinates of the corresponding variants.
It requires +/- 20 minutes to iterate through the entire dbSNP VCF file, which is by far the longest step of the build procedure.
Since ALFRED is at its end-of-life and will not be updated, I'm suspended best practice and dropping the computed microhap coordinates here as a "curated" input, rather than re-computing these values from scratch each time as part of the build.
