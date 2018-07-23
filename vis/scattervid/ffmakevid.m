## ffmakevid.m
##
## Makes a video from a framecache (collection of consecutively named
## stillframe images). Uses ffmpeg to encode the video.  ffmpeg must be
## installed on the system.  On Ubuntu this can be installed with
## 'sudo apt-get install ffmpeg'.
##
## We set the options to produce an h.264 encoded .mp4 file that will
## be playable by Apple Quicktime players.  (Future modifications may
## facilitate selecting other formats, but this suits my needs.)
##
function ffmakevid(outfile,         # Outfut file name, e.g. "movie.mp4"
                   framecachepat,   # Framecache pattern, e.g. "framecache/FIG__%04d.png" 
                   varargin         # Other args as "Prop",Value pairs
                  )                 #  (Warning: unrecognized or misspelled props
                                    #   will be silently ignored.)

  pixwidth  = parsevararg("pixwidth", 1000, varargin{:})
  framerate = parsevararg("framerate", 10, varargin{:})
  bitrate   = parsevararg("bitrate", 2000000, varargin{:})

  syscmd = cstrcat("ffmpeg",
            " -r ", num2str(framerate),   # Input frame rate (otherwise assumes)
            " -i ", framecachepat,        # Input file(s)
            " -f mp4",                    # Mp4 container, h264 implied I think
            " -c:v libx264 ",             # h.264 encoding
            " -pix_fmt yuv420p",          # Quicktime gets confused if 422.
            " -preset veryslow ",         # Slow encode but better result (x264)
            " -tune animation ",          # Tweaks for animation input (x264)
            " -r ", num2str(framerate),   # Output frame rate
            " -b:v ", num2str(bitrate,"%d"),# Bits per second 
            " -vf 'scale=", num2str(pixwidth), ":trunc(ow/a/2)*2'",
            " ", outfile                  # Destination file
           );

  syscmd # Output command-line to stdout

  system(syscmd);

end


##
function ret = parsevararg(fieldname, defvalue, varargin)
  found=0;
  for i=1:length(varargin);
    if (strcmp(toupper(varargin{i}),toupper(fieldname)))
      found=i;
      break;
    end
  end
  if (found>0)
    ret = varargin{found+1};
  else
    ret = defvalue;
  end
end
