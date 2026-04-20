# ns-3 THz-NTN Module

> **World's First ns-3 Module for Terahertz Non-Terrestrial Networks with Deep Satellite + mmWave + O-RAN Integration**

[![ns-3](https://img.shields.io/badge/ns--3-3.43-blue)](https://www.nsnam.org/)
[![License: GPL-2.0](https://img.shields.io/badge/License-GPL--2.0-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.html)
[![Tests](https://img.shields.io/badge/tests-12%2F12%20passing-green)]()

Research-grade extension to the `ns3-ntn-toolkit` that adds full **Terahertz (100 GHz -- 10 THz)** capabilities to the integrated satellite/mmWave/O-RAN simulation stack. Designed in collaboration with UNLab (Northeastern University) for use with sub-THz missions such as **TeraLink-1** (NSF Award #2346487).

---

## Key Features

| Capability | Description |
|---|---|
| **HITRAN molecular absorption** | Altitude-stratified atmosphere (ITU-R P.835) with H2O and O2 absorption lines 100 GHz -- 10 THz |
| **Weather attenuation** | ITU-R P.838/P.840 rain, fog, snow, dust models extended to THz |
| **Scintillation** | ITU-R P.618 amplitude + phase scintillation with AR(1) time series |
| **Pointing error** | Vibration, J2 perturbation, atmospheric refraction, tracking latency |
| **Hardware impairments** | Rapp / Saleh PA, Lorentzian phase noise, ADC SQNR, I/Q imbalance, DPD |
| **Ultra-Massive MIMO** | UPA / UCA / Cassegrain arrays up to 128x128 elements with beam squint |
| **Hierarchical beamforming** | DFT codebook, multi-resolution search, hybrid analog + digital |
| **EKF beam tracking** | Satellite ephemeris-assisted Kalman filter for LEO pointing |
| **Inter-Satellite Links** | Vacuum propagation at THz with Dijkstra routing and Doppler |
| **Reconfigurable Intelligent Surfaces** | Space-borne / aerial / ground RIS with N^2 SNR scaling |
| **ISAC** | Joint radar-communication for space debris (cm-level resolution) |
| **Novel waveforms** | OFDM / DFT-s-OFDM / OTFS / AFDM / SC-FDE for high-Doppler LEO |
| **Atmospheric windows** | 140 / 220 / 340 / 410 / 460 GHz sub-THz windows |
| **O-RAN integration** | THz-aware xApps (beam management, spectrum, RIS, ISAC) |

## Deep Module Integration (Not Standalone)

Unlike hobbyist THz modules, every value in `thz-ntn` flows through the **actual satellite and mmWave module APIs**:

- `ThzNtnFreeSpaceLoss` **extends `SatFreeSpaceLoss`** -- drop-in replacement for the satellite channel pipeline
- `ThzNtnPhy` uses actual `SpectrumValue` PSDs through `mmWaveInterference`
- `ThzNtnPhySat` delegates MCS selection to `SatWaveformConf::GetBestWaveformId()`
- `ThzNtnMac` filters sub-bands using `ThzNtnMolecularAbsorption` with actual `MobilityModel` positions
- `ThzNtnAntennaArray` can delegate to `SatAntennaGainPattern` for beam patterns
- `ThzNtnBeamTracking` EKF state updated from satellite `MobilityModel::GetVelocity()`
- Doppler / distance / elevation / altitude all derived from real `MobilityModel` 3D positions

**No hardcoded values.** Every loss, gain, capacity, and metric is computed from actual physics using actual parameters.

---

## Repository Layout

```
contrib/thz-ntn/
├── model/                      # 48 source files (24 .h + 24 .cc)
│   ├── thz-ntn-channel-model.{h,cc}
│   ├── thz-ntn-free-space-loss.{h,cc}    # Extends SatFreeSpaceLoss
│   ├── thz-ntn-molecular-absorption.{h,cc}
│   ├── thz-ntn-weather-attenuation.{h,cc}
│   ├── thz-ntn-scintillation.{h,cc}
│   ├── thz-ntn-pointing-error.{h,cc}
│   ├── thz-ntn-hardware-impairments.{h,cc}
│   ├── thz-ntn-spectrum.{h,cc}
│   ├── thz-ntn-link-budget.{h,cc}
│   ├── thz-ntn-phy*.{h,cc}                # Base / Sat / Ground
│   ├── thz-ntn-waveform.{h,cc}            # OFDM / OTFS / AFDM / ...
│   ├── thz-ntn-antenna-array.{h,cc}       # UM-MIMO / Cassegrain
│   ├── thz-ntn-beamforming.{h,cc}         # Hierarchical DFT codebook
│   ├── thz-ntn-beam-tracking.{h,cc}       # EKF predictor
│   ├── thz-ntn-mac*.{h,cc}
│   ├── thz-ntn-isl-{channel,link}.{h,cc}
│   ├── thz-ntn-ris*.{h,cc}
│   └── thz-ntn-isac*.{h,cc}
├── helper/                     # 2 source files
│   └── thz-ntn-helper.{h,cc}
├── examples/                   # 8 scenario scripts (in-tree) +
│                                # scratch/thz-ntn-demo.cc (runnable)
├── test/                       # 1 test suite with 12 tests
├── doc/                        # Detailed docs (this directory)
└── CMakeLists.txt
```

---

## Quick Start

### 1. Build

```bash
cd ns-3-dev
./ns3 configure --enable-modules=thz-ntn
./ns3 build thz-ntn
```

### 2. Run tests

```bash
./ns3 run "test-runner --suite=thz-ntn --verbose"
```

Expected: `PASS thz-ntn 0.030 s` with 12 passing tests.

### 3. Run the demo

Copy `examples/thz-ntn-demo.cc` to `scratch/` (or use the ready-made scratch version), then:

```bash
./ns3 build scratch/thz-ntn-demo
./build/scratch/ns3.43-thz-ntn-demo-debug                 # All 8 scenarios
./build/scratch/ns3.43-thz-ntn-demo-debug --example=1      # Just scenario 1
```

CSV datasets land in `thz-ntn-results/` (see [doc/EXAMPLES.md](doc/EXAMPLES.md) for details).

---

## Documentation

- [doc/EXAMPLES.md](doc/EXAMPLES.md) -- Walk-through of all 8 demo scenarios with expected values
- [doc/MODULE_REFERENCE.md](doc/MODULE_REFERENCE.md) -- Per-class reference (APIs, attributes, formulas)
- [doc/INTEGRATION.md](doc/INTEGRATION.md) -- How this module hooks into the satellite and mmWave pipelines
- [doc/VALIDATION.md](doc/VALIDATION.md) -- Verified against analytical formulas and ITU-R models

---

## Minimal Usage Example

```cpp
#include <ns3/thz-ntn-link-budget.h>
#include <ns3/thz-ntn-free-space-loss.h>
#include <ns3/thz-ntn-molecular-absorption.h>
#include <ns3/constant-position-mobility-model.h>

using namespace ns3;

// Set up the full loss stack
auto absorption = CreateObject<ThzNtnMolecularAbsorption>();
auto fsl        = CreateObject<ThzNtnFreeSpaceLoss>();
fsl->SetMolecularAbsorptionModel(absorption);
fsl->EnableMolecularAbsorption(true);

// Link budget with actual positions
auto ground = CreateObject<ConstantPositionMobilityModel>();
ground->SetPosition(Vector(6371000, 0, 0));
auto sat    = CreateObject<ConstantPositionMobilityModel>();
sat->SetPosition(Vector(6371000 + 550000, 0, 0));

auto lb = CreateObject<ThzNtnLinkBudget>();
lb->SetFreeSpaceLossModel(fsl);
auto result = lb->ComputeTeraLinkBudget(sat, ground);

std::cout << "FSPL:         " << result.fspl_dB              << " dB\n";
std::cout << "Absorption:   " << result.molecularAbsorption_dB<< " dB\n";
std::cout << "SNR:          " << result.snr_dB               << " dB\n";
std::cout << "Capacity:     " << result.shannonCapacity_Gbps << " Gbps\n";
```

At 225 GHz, LEO 550 km, zenith: **FSPL=194.30 dB** (matches `20*log10(4*pi*d*f/c)` exactly).

---

## Tested Environment

- ns-3.43 on Ubuntu 24.04
- GCC 14, CMake 3.27
- Eigen3, Boost, GSL, SQLite3

## Citation

```bibtex
@misc{uzair2026thzntn,
  author       = {Muhammad Uzair},
  title        = {ns-3 THz-NTN Module: Terahertz Non-Terrestrial Networks
                  with Integrated Satellite and mmWave Simulation},
  year         = {2026},
  howpublished = {\url{https://github.com/Muhammaduazir69/ns3-thz-ntn}}
}
```

## License

GPL-2.0-only (same as the satellite module it extends).

## Acknowledgements

- UNLab at Northeastern University (Prof. J. M. Jornet) for `TeraSim` inspiration and THz-NTN research direction
- Magister Solutions for SNS3 satellite module
- NYU Wireless for the ns-3 mmWave module
- ITU-R, 3GPP TR 38.811, and ESA D-band satellite links project for propagation models

---

**Muhammad Uzair** - 2026 - Independent Researcher - 5G/6G NTN
