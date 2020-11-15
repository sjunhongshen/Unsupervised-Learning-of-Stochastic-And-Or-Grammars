# Sketch Face Testcase

## 1.Usage

Executable arguments:
- `-n` (required): specify number of samples
- `-g` (required): specify the number of gibbs sampling
- `-o` (optional): specify the output directory, default is `face_[Num of gibbs]`
- `-b` (optional): if used, then samping negative samples

```bash
# output directry is faces_200/, and sampling 100 positive sample
./cartoon_face -n 100 -g 200
# output directry is neg-samples/, and sampling 50 negative sample
./cartoon_face -n 50 -o neg-samples -b
```

