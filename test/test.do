version 14.1
set seed 5371
set linesize 112

cap which multe
if _rc {
    disp as err "Please install multe"
    exit 111
}

capture program drop main
program main
    qui do test/test_replicate.do
    qui do test/test_unit.do
    qui do test/test_weights.do

    multe_unit_tests
    multe_replicate_tests
    multe_weight_tests
end

main
