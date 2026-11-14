# Bare finite element reference inputs

Credit: Gong Xiangyu, He Yanli and Li Xiongyan, AirMem response-side screening, numerical inputs carried forward from the V18.1 reproducibility archive. These numerical model inputs are provided under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/), as identified in the source archive. This distribution retains the original files.

- `full_model_pure.mac`: model mesh geometry, membrane and cable connectivity.
- `FE_modal_sensor_z585.csv`: bare model frequencies and vertical components at the 12 sensors.
- `my_freqs_585z.npy`: bare modal frequencies used in the frequency library.
- `modes600.npz`: mass-normalized full-field displacement vectors from the bare prestressed model.

The files are model inputs, separate from the public monitoring records. They do not contain measured screening results, fitted added mass, or corrected frequencies.
