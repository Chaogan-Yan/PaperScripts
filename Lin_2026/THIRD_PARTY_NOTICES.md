# Third-party software

This repository does not redistribute DPABI, DPABISurf, ComBat, MATLAB, or
their template and atlas files. The MATLAB scripts call functions supplied by
locally installed dependencies. Users must obtain those dependencies
separately and comply with their licenses.

- [DPABI/DPABISurf](https://github.com/Chaogan-Yan/DPABI) is distributed under
  the GNU Lesser General Public License v2.1. The FC reporting script uses
  DPABI I/O, reslicing, atlas, and GRF-compatible functionality.
- [ComBatHarmonization](https://github.com/Jfortin1/ComBatHarmonization)
  provides the MATLAB ComBat implementation called by
  `code/combat_falff_maps.m`; its MATLAB code is distributed under the MIT
  License.
- MATLAB and the named MathWorks toolboxes are proprietary dependencies and
  are not included here.

The DIRECT imaging data, participant-specific labels, DPABI templates, atlases,
and all study-derived participant-level files are also excluded from this
repository.
