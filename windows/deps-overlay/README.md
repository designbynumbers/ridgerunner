# deps-overlay

CMake + OpenMP/TLS patches for sibling builds when CI checks out stock
designbynumbers tags. Local `../ridge_tsnnls` and `../ridge_plcurve` already
contain these changes on `feature/Build_Windows_Multi_Thread`.

Apply after checkout:

```bat
xcopy /E /Y windows\deps-overlay\tsnnls\* ..\ridge_tsnnls\
xcopy /E /Y windows\deps-overlay\plcurve\* ..\ridge_plcurve\
```

(CI does the equivalent before `run_all.py`.)
