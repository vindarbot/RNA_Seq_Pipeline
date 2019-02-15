#!/bin/bash

### VARIABLES ###

EXP="Experience"

SAMPLES=$(ls "$EXP")

NB_SAMPLES=$(ls "$EXP"/*/*R1* | wc -l)

