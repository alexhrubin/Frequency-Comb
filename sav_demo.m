%% Yooo this is gonna walk thru the changes to Waveguide.m and show how
%% to use DualWaveguide.m. I'll break it up into sections,
%% I'll describe each section, then you should run it to see what it does.

%% Remember we were having issues before, because it was annoying
%% to update the shape of the waveguide in Waveguide.m?
%% We had the functions describing the waveguide edge hard-coded inside
%% the class, and when we updated them we would forget to update the\
%% derivatives, which also had to be done by hand. Now, when we define
%% a waveguide, we explicitly give it two function handles, one for the
%% bottom edge, one for the top edge. We don't have to think about the
%% derivatives at all, matlab calculates them for us automagically.
%% Here I will make 2 function handles: one for the lower edge, one
%% for the upper. Then I'll make a waveguide from them and do some basic
%% stuff that you've seen before, just to show it works.

upper_edge = @(z) 0.075 * exp(0.025 * z);
lower_edge = @(z) -0.075 - 0.002 * z;

wg = Waveguide(12, 2*pi/1.55, lower_edge, upper_edge);

wg.visualize_waveguide([0 40])
wg.plot_eigenmodes('all', 0)
wg.plot_eigenmodes('all', 40)

%% OK. Everything else about the Waveguide class is basically the same.
%% Now the DualWaveguide class. We define this exactly the same way as
%% the Waveguide class. We create 4 function handles, one for each edge
%% of the two waveguides that we have. Once we do that, we can visualize
%% the two waveguides in the same way as we usually do. I will refer to 
%% the lower waveguide as wg1, and the upper one as wg2.
%% Now, we have to pay attention to two things.
%% (1) For the math we did to work, we are assuming that the bottom edge
%% of the bottom waveguide stays right on top of the z-axis and does not
%% change.
%% (2) When we create the DualWaveguide, we specify the distance between
%% the centers of the waveguides at z=0.
%% I'll make an example and when it's visualized it will be clearer.

clf

wg1_lower_edge = @(z) 0;   %% the bottom edge of the bottom waveguide stays 0.
wg1_upper_edge = @(z) 0.15 + 0.002 * z;

wg2_lower_edge = @(z) -0.002 * z;
wg2_upper_edge = @(z) 0.15 + 0.002*z;

sep = 0.5;

dwg = DualWaveguide(12, 2*pi/1.55, sep, wg1_lower_edge, wg1_upper_edge, wg2_lower_edge, wg2_upper_edge);
dwg.visualize_waveguide([0 40])

%% So we have two waveguides each expanding. Now we have some functions
%% that let us solve for the betas, and we can also plot the eigenmodes.
%% You will see that at z=0, the waveguides are the same size, so the 
%% mode looks like two identical humps,one in each waveguide.

dwg.getbeta_position(0)
dwg.plot_eigenmodes(0, 'all')

%% If we now move further away to z=40, the waveguides are not the same
%% width anymore, so the mode does something more complicated.

dwg.plot_eigenmodes(40, 'all')

%% Unfortunately there are some difficulties with solving for the betas,
%% you can see that here: Two of our betas here are identical, which
%% obviously should not happen... It's easy to solve for the betas
%% in a specific case, but it is tricky to come up with a method which
%% guarantees you will always find all of them, and never duplicate them.
%% That is something you can try to work on improving.

dwg.getbeta_position(40)

%% Here is why that is tricky: what is happening we solve for the betas,
%% is we are finding the zeros of the function dual_wg_eigenproblem.m
%% (there is a copy of that function in DualWaveguide.m)
%% dual_wg_eigenproblem.m is a really weird wavy function that changes
%% shape a lot, so the zeros are hard to predict. One solution, suggested
%% by Weiliang, is that the betas of the two waveguide system will be
%% close to the betas of the single-waveguide, so we can use those betas
%% as a starting point and let matlab search nearby. Unfortunately, sometimes
%% there are two betas close to each other, so making sure you hit both is tough.
%% Let's look at the plot for the dual waveguide we just made, as an example.

a1 = 2 * dwg.wg1.half_width(40);
a2 = 2 * dwg.wg2.half_width(40);
sep = dwg.edge_sep(40);

f = @(b) dual_wg_eigenproblem(b, 2*pi/1.55, 12, a1, a2, sep);

beta_range = 0:0.01:14;

clf
plot(beta_range, arrayfun(f, beta_range));
grid on

%% See there are 3 distinct zeros? We got three betas dwg.getbeta_position(40),
%% but two are identical, meaning we missed the third one. This code
%% DOES work most of the time, but this is a case showing it's not perfect...

%% Anyway, once we get the betas, we run dwg.boundary_conds, which
%% solves for the boundary conditions (this part I know is correct),
%% and then we have the full eigenmode function to do whatever we want with.

%% You might wanna mess around with the DualWaveguide class, make some
%% weird shapes and see what happens. Try to break stuff and see if you can
%% fix it.

%% Also look at dual_wg_overlap.m. It is the extension of overlap.m
%% to work for two waveguides (obviously). The overlap integral
%% really only depends on where the edges of the waveguides are, so
%% now that there are four edges instead of two, the function had to be
%% changed. I am pretty sure that I did it correctly, but I have not
%% tried to make a full simulation with ode45 yet! Of course that is
%% where we ultimately need to go.

