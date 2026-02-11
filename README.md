JFEM is a linear static finite element solver written in Julia that reads MSC Nastran BDF (Bulk Data File) input decks and produces results compatible with Nastran SOL 101. It is designed for educational use and as a validation tool for structural analysis workflows.

Key features:

Reads standard fixed-format and free-format Nastran BDF files
Solves linear static analysis (equivalent to SOL 101)
Exports results as VTK files (for ParaView visualization) and JSON
Includes a validation tool for comparing results against Nastran F06 output
