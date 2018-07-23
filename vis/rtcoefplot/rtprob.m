function rtprob(fname)

  DATA = load(fname);

  sini = DATA(:,1);
  PROBS = DATA(:,2:7);
  TPROBS = sum(PROBS,2);
  NPROBS = PROBS ./ TPROBS;  # Normalized probabilities

  RP  = NPROBS(:,1);
  TP  = NPROBS(:,2);
  RSV  = NPROBS(:,3);
  TSV  = NPROBS(:,4);
  RSH  = NPROBS(:,5);
  TSH  = NPROBS(:,6);

  XX = asin(sini);

  clf();
  hold("on");
  plot (XX, RP, "r");
  set(get(gca()).children(end), "linewidth", 2.0);
  plot (XX, TP, "m");
  set(get(gca()).children(end), "linewidth", 0.5);
  plot (XX, RSV, "b");
  set(get(gca()).children(end), "linewidth", 2.0);
  plot (XX, TSV, "c");
  set(get(gca()).children(end), "linewidth", 0.5);
  plot (XX, RSH, "g");
  set(get(gca()).children(end), "linewidth", 2.0);
  plot (XX, TSH, "y");
  set(get(gca()).children(end), "linewidth", 0.5);


end
