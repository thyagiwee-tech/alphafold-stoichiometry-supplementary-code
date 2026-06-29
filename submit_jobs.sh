find . -maxdepth 1 -type d \( ! -name . \) -exec bash -c "cd {} && qsub af1.sh" \;
~                                                                                        