# Drop40 nine-frame data correction

The final Drop40 temperature series is:

| Frame | Full scan | Nominal temperature (°C) | Measured `temp3` (°C) |
|---:|---:|---:|---:|
| 193 | 587193 | 40 | 39.9 |
| 194 | 587194 | 50 | 49.9 |
| 195 | 587195 | 55 | 54.8 |
| 196 | 587196 | 60 | 59.9 |
| 197 | 587197 | 65 | 64.5 |
| 198 | 587198 | 70 | 69.6 |
| 199 | 587199 | 75 | 74.5 |
| 200 | 587200 | 80 | 79.2 |
| 201 | 587201 | 150 | 145.4 |

An earlier eight-frame export, kept only in the private working tree, contained
frames 193-200 and its fitting configuration incorrectly associated frame 199
with 80 °C and frame 200 with 150 °C. Every result produced from it is
superseded.

The corrected inputs are committed in `example_data/drop40/`:

- `drop40_FR_9frames_norm.txt`
- `drop40_IP_9frames_norm.txt`
- `drop40_OOP_9frames_norm.txt`
- `drop40_9frame_manifest.csv`

They can be regenerated from the original sector CSV exports, which are not in
this repository, with:

```bash
python scripts/preparation/build_drop40_nine_frame_inputs.py \
  --data-root /path/to/external_data \
  --output-dir example_data/drop40
```

The source exposure time and transmission are constant across the nine scans.
The norm-corrected files multiply each scan by
`median(norm across all scans) / scan norm`. This preserves approximately the
original intensity scale while correcting the roughly 6% recorded monitor
variation. Raw and normalized results should be compared before interpreting
areas quantitatively.

All peak-fitting outputs, spreadsheets, and slides created before this
correction are superseded and must not be used as final numerical results.
