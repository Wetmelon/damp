# Motor FOC (opt-in)

FOC current regulation and carrier PWM for three-phase VSI. Not pulled by
`control.hpp` — include what you need:

```cpp
#include "damp/motor/foc.hpp"         // FOController → clamped Vdq
#include "damp/motor/modulation.hpp"  // SVPWM / SPWM / DPWM → duties
#include "damp/motor/spm.hpp"         // surface-PM Kt/λ helpers (optional)
```

Teaching demo: `examples/motor/foc/` (nameplate → FOC tick → host SIL plots).

This is the FOC slice. Drive / machine adapters / full motor pack land later.
