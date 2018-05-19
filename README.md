__!! WARNING : WORK IN PROGRESS !!__

# nanoes
Nano EcmaScript, easily embeddable single header C++ 14 runtime.

## Licence
MIT

## Roadmap

### Status

Done:
+ lexer
+ parser / script-gen
+ GC
+ some basic arithmetic operators
+ function (s)
+ if statements
+ while loops

TODO:
+ objects
+ arrays
+ full operator set
+ for loops
+ var variables

### Language

This runtime does not aim for perfect per-bug compatibility but put embeddability and sane execution first. (See the design goals)
With the exception of additions involving strings the runtime may rather fault than try to do operations on incompatible types with unpredictable results.


## Design

Design priorities in order of importance.
- Embeddability first (easy access of generic C++ values)
- Easy to build (single header+sourcefile to add)
- Small code (Compact maintainable codebase)
- Low/configurable memory utilization.
- Decent performance if possible (not a primary target, focus on game-like workloads)
