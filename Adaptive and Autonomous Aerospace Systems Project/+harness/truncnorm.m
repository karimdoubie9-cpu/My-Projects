function x = truncnorm(mu, sigma, range, sz)
%TRUNCNORM  Normal draws resampled until they fall inside [range(1) range(2)].
%
%   x = harness.truncnorm(0.4, 0.08, [0.15 0.65], [n 1])
%
%   Resampling, not clipping. Clipping would pile probability mass onto the two
%   bounds and put spurious spikes at the edges of every histogram; resampling
%   leaves the interior shape of the normal intact.
%
%   The bounds are physical rather than cosmetic: a package cannot have negative
%   mass, and a rotor cannot deliver MORE thrust than it was commanded, which is
%   why the effectiveness range is hard-capped at 1.
%
%   SHARED ON PURPOSE. harness.montecarlo and run_montecarlo both call this, so
%   both consume the random stream in exactly the same order and produce
%   identical draws from the same seed. Keeping a private copy in either place
%   would let the two campaigns silently diverge the first time a distribution
%   was retuned - and the whole point of the second runner is that its numbers
%   match the first.
%
%   Do not reorder or add randn calls here without accepting that every stored
%   campaign becomes irreproducible.

arguments
    mu    (1,1) double
    sigma (1,1) double
    range (1,2) double
    sz    (1,2) double
end

x     = mu + sigma*randn(sz);
bad   = x < range(1) | x > range(2);
guard = 0;

while any(bad(:)) && guard < 1000
    x(bad) = mu + sigma*randn(nnz(bad), 1);
    bad    = x < range(1) | x > range(2);
    guard  = guard + 1;
end

x = min(max(x, range(1)), range(2));   % belt and braces after the guard limit
end
