# combine.m
#
# Combines output from two seis files, presumed to represent different
# synthetic runs for the same seismometer in the same model and sim
# conditions, to produce an aggregate seis file with a higher phonon
# count.
#
# Also works to combine out_mparams.octv file to produce one
# representing the correct summed phonon count.
#
function combine(sfn1, # seis filename 1
                 sfn2, # seis filename 2
                 osfn  # output seis filename
                )

  S1 = load(sfn1);
  if (isfield(S1,"SEIS")) S1=S1.SEIS; end # Catch nested format
  S2 = load(sfn2);
  if (isfield(S2,"SEIS")) S2=S2.SEIS; end # Catch nested format

  SEIS = S1;  # This will be the output

  # TODO: Check for compatibility between Location parameters, etc.
  # to confirm that these files really should be combined

  if (isfield(S1,"TraceXYZ"))
    SEIS.TraceXYZ = S1.TraceXYZ + S2.TraceXYZ;
    SEIS.TracePS = S1.TracePS + S2.TracePS;
    SEIS.CountPS = S1.CountPS + S2.CountPS;
  end
  if (isfield(S1,"NumPhonons"))
    SEIS.NumPhonons = S1.NumPhonons + S2.NumPhonons;
  end

  save("-text", osfn, "SEIS");

end
