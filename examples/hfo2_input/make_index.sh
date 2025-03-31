#!/bin/bash
echo -e "a" {65..96} "\nname 3 bottom\na" {1..64} {97..165} "\nname 4 top\nq" | gmx make_ndx -f hfo2.gro