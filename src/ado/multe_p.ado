capture program drop multe_p
program multe_p, sclass sortpreserve
    version 14.1
    syntax [anything] [if] [in], [*]
    disp as err "{bf:error:} postestimation not yet available"
    error
end
