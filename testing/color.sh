#!/bin/bash

# Correct usage of echo with escape sequences
printf "\e[31mThis is red text\e[0m\n"
echo -e "\e[1;32mThis is green bold text\e[0m"
echo -e "\e[4;33mThis is underlined yellow text\e[0m"

