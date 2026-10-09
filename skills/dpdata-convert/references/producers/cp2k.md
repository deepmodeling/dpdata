# CP2K producer contract

Use this reference when converting CP2K `*.out`, AIMD, or trajectory output.
Select the dpdata format that matches the producer output (`cp2k/output` or
`cp2k/aimd_output`) and record the CP2K release and the output sections that
were enabled.

Before accepting the converted data, verify:

- the coordinate and cell sections describe the same frames as the energy and
  force sections;
- the calculation's periodicity and cell vectors are preserved;
- energies, forces, and any stress/virial values use the units and sign
  convention expected by the downstream consumer; and
- restarts or truncated runs did not create unmatched frames.

CP2K output details vary by release and input print settings. Treat the CP2K
manual and the output producer version as the authority. Parser defects belong
in the parser project; this reference defines the validation questions to ask.
