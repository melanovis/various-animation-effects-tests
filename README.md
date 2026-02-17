# various-animation-effects-tests
various test scripts for some animation effects I thought were cool and decided to try out

- Portal 2 latex text effect: A very simple script which splits an image up into a series of tiles and then rotates them around using a transformation matrix. The rotations are set up to cascade diagonally with a tuned second order response to give them a natural looking transition with a tiny bit of overshoot. This was originally developed as a way to transition between latex formulas for a uni project but was scrapped.
- Markov chain annealing text effect: This short script takes a series of cells of a given index and tries to markov-chain anneal them into a specified map (in this case an image of a word where each letter corresponds to some index). This has been (mostly) successful and the transition effect is unique and elegant so i'd like to use it in a future project.

No LLM-derived code is used in this work.
