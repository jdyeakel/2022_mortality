using(DataFrames)
using(CSV)
using(RCall)
using(LinearAlgebra)
using(Distributions)
using(Distributed)
using(UnicodePlots)

haywardfulldata = CSV.read("[YOUR FILE PATH]/data/data_hayward_all.csv",header=true,DataFrame);

# "Panthera leo"
# "Crocuta crocuta"
# "Panthera pardus"
# "Cuon alpinus"
# "Lycaon pictus"
# "Acinonyx jubatus"
# "Panthera tigris"

#Caculcate mean preferred mass per predator
preds = unique(haywardfulldata[!,:Predator]);

function genpredprey()
preyinds = Array{Float64}(undef,0);
predinds = Array{Float64}(undef,0);
for i=1:length(preds)
    pref = haywardfulldata[!,:PercentOfKills][findall(x->x==preds[i],haywardfulldata[!,:Predator])];
    pref[findall(isnan,pref)].=0.;
    pref = round.((pref ./ sum(pref))*1000); #10000
    preyind = sum(pref);

    meanpredmass = mean(haywardfulldata[!,:Predbodymasskg][findall(x->x==preds[i],haywardfulldata[!,:Predator])]);
    predmassSD = 0.25*meanpredmass;
    predbodysizedist = Normal(meanpredmass,predmassSD);

    #DRAW BODY SIZES FOR PREDATOR INDIVIDUALS
    predinds_draw = abs.(rand(predbodysizedist,Int64(preyind)));
    predinds = [predinds; predinds_draw];
    preyi = haywardfulldata[!,:Prey][findall(x->x==preds[i],haywardfulldata[!,:Predator])];
    
    for j=1:length(preyi)
        if pref[j] > 0.
            #draw body masses
            #Preybodymasskg
            #Preybodymasskg34adultfemalemass
            meanmass = haywardfulldata[!,:Preybodymasskg34adultfemalemass][findall(x->x==preds[i],haywardfulldata[!,:Predator])][j];
            massSD = 0.25*meanmass;
            preybodysizedist = Normal(meanmass,massSD);
            numbers = Int64(pref[j]);
            preyinds_draw = abs.(rand(preybodysizedist,numbers));
            preyinds = [preyinds; preyinds_draw];
        end
    end
end
return(predinds,preyinds)
end



function meanpredpreysize(predinds,preyinds,sizebins)
    maxprey = log10(maximum(preyinds));
    minprey = log10(100); #Run this only for prey = 100 KG to max KG
    stepsize = (maxprey-minprey)/sizebins;
    preysizeclass = collect(minprey:stepsize:maxprey);
    meanpredsize = Array{Float64}(undef,(length(preysizeclass)-1));
    for i=2:length(preysizeclass)
        preysizeinds_pos = findall(x->((x>10^preysizeclass[i-1]) && (x < 10^preysizeclass[i])), preyinds);
        if length(preysizeinds_pos) > 0
            meanpredsize[i-1] = mean(predinds[preysizeinds_pos]);
        else 
            meanpredsize[i-1] =NaN;
        end
    end
    filledspots = findall(!isnan,meanpredsize)

    return meanpredsize[filledspots], (10 .^preysizeclass[filledspots])
end


#Do this a bunch of times to get a number of replicates
sizebinsvec = collect(90:100);

sizebinsvec = repeat([100],inner=1000)
reps = length(sizebinsvec);


predpreyreps = Array{Float64}(undef,reps,maximum(sizebinsvec),2);
for r = 1:reps
    sizebins = sizebinsvec[r];
    predinds,preyinds = genpredprey();
    meanpredmass,meanpreymass = meanpredpreysize(predinds,preyinds,sizebins);
    lsp = size(meanpredmass)[1];
    predpreyreps[r,1:lsp,1] = meanpreymass;
    predpreyreps[r,1:lsp,2] = meanpredmass;

    predpreyreps[r,(lsp+1):maximum(sizebinsvec),1] .= NaN;
    predpreyreps[r,(lsp+1):maximum(sizebinsvec),2] .= NaN;
end

preyvec = vec(predpreyreps[:,:,1]);
preyvalues = preyvec[findall(!isnan,preyvec)];
predvec =  vec(predpreyreps[:,:,2]);
predvalues = predvec[findall(!isnan,predvec)];

scatterplot(log.(preyvalues),log.(predvalues))

R"""
model5 = lm(log($(predvalues)) ~ log($(preyvalues)))
summary(model5)
fitintercept = model5[[1]][[1]];
fitslope = model5[[1]][[2]];
CI = confint(model5,level=0.95)
intlow = CI[[1]];
slopelow = CI[[2]];
inthigh = CI[[3]];
slopehigh = CI[[4]];
"""
fitintercept = @rget fitintercept;
fitslope = @rget fitslope;
fitinterceptlow = @rget intlow;
fitslopelow = @rget slopelow;
fitintercepthigh = @rget inthigh;
fitslopehigh = @rget slopehigh;


#Export for the mathematica notebook
fit_table = DataFrame([fitintercept fitinterceptlow fitintercepthigh; fitslope fitslopelow fitslopehigh], :auto);
rename!(fit_table,[:Fit,:FitLow,:FitHigh])

CSV.write("[YOUR FILE PATH]/data/ppmr_fit_table_revreps.csv",fit_table; header=true);

predpreysizetable = DataFrame([preyvalues predvalues],:auto);
rename!(predpreysizetable,[:preymass,:predmass])

CSV.write("[YOUR FILE PATH HERE]/data/predpreymass_tablereps.csv",predpreysizetable; header=false);
