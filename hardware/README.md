# hardware/

## MIMO_LCL_conditioning.net (Supplementary Data S11)

Netlist of the six-channel conditioning model used for the circuit analysis of Section 2.1.2 of the paper. It was
exported from LTspice in ExpressPCB netlist format (part IDs, net names and net connections tables).

**What it contains**

- Source side: three 325 V peak, 50 Hz sources at 0°, 120° and 240° (`Vinv1`, `Vbinv1`, `Vcinv1`) behind an
  LCL output filter (`L1`–`L3` 2 mH, `Cfa1`–`Cfc1` 10 µF with `Rd1`–`Rd3` 10 Ω damping, `L4`–`L6` 1 mH), and the
  ideal grid-side sources `Va1`, `Vb1`, `Vc1`.
- Per voltage channel (A, B, C): potential-transformer divider `R_pt` 65 kΩ / 1 kΩ; divider `R1` 10 kΩ over `R2`
  3.3 kΩ; coupling capacitor `C1` 1 µF; bias resistors 1 MΩ / 1 MΩ to 3.3 V and ground; `Cbias` 100 nF
  (`C3`, `C6`, `C9`); anti-aliasing `Rf` 1 kΩ with 100 nF (`C2`, `C5`, `C8`); ADC input load `R_load` 100 MΩ.
- Per current-proxy channel (A, B, C): `R_series` 3.3 kΩ, `Cproxy` 1 µF, and 1 MΩ / 1 MΩ bias.
- Supply: `V1` 3.3 V with `C_bulk1` 10 µF and `C_bypass1` 10 nF.

**What it does not contain**

- The operational-amplifier symbols. The analysis used Microchip's MCP6001 SPICE macromodel (the single-channel
  device electrically equivalent to one half of the fitted MCP6002), which is Microchip's and is not
  redistributed here; it is available from Microchip's MCP6001 product page. The op-amp nets are still named in
  the netlist (`V_signalX` at each first-buffer input, `OpAmp_OutX` at its output, `V_to_ADC_X` and `I_to_ADC_X`
  at the converter inputs), so each unity-gain follower is re-inserted between those nets.
- The as-built wiring faults. The netlist models all three voltage chains complete, as designed. On the bench
  board the Va and Vc first-buffer inputs are unconnected (Section 2.5 of the paper).

**Edit made for release.** The exported file encoded the micro prefix as a byte that did not survive transfer; it
has been replaced with the SPICE prefix `u` (so `1u` is 1 µF). No value or connection was changed.

The analytic derivation of the chain gain, phase and poles from these values is Supplementary Document S1.
