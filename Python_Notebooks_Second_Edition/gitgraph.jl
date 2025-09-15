import Downloads
import ZipFile

using SHA
using CSV
using DataFrames
using Graphs
using GLM
using JSON3
using Statistics

git_zip = "git_web_ml.zip"
if !isfile(git_zip)
    Downloads.download("https://snap.stanford.edu/data/git_web_ml.zip", git_zip)
end

function ingest_to_df(archive::ZipFile.Reader, filename::AbstractString)
    idx = only(findall(x -> x.name == filename, archive.files))
    return CSV.read(read(archive.files[idx]), DataFrame)
end

git_archive = ZipFile.Reader(git_zip)
edges_df = ingest_to_df(git_archive, "git_web_ml/musae_git_edges.csv");
classes_df = ingest_to_df(git_archive, "git_web_ml/musae_git_target.csv");
idx = only(findall(x -> x.name == "git_web_ml/musae_git_features.json", git_archive.files))
attributes = JSON3.read(read(git_archive.files[idx]))
close(git_archive)

edges_df .+= 1
classes_df.id .+= 1

# I picked feature 1793 as protected, but maybe we will like some other
classes_df.protected = [1793 in last(x) for x in attributes]

gh = SimpleGraph(nrow(classes_df))
for (src, dst) in eachrow(edges_df)
    add_edge!(gh, src, dst)
end

function deg_class(gh, class)
    deg_ml = zeros(Int, length(class))
    deg_web = zeros(Int, length(class))
    for edge in edges(gh)
        a, b = edge.src, edge.dst
        if class[b] == 1
            deg_ml[a] += 1
        else
            deg_web[a] += 1
        end
        if class[a] == 1
            deg_ml[b] += 1
        else
            deg_web[b] += 1
        end
    end
    return (deg_ml, deg_web)
end

classes_df.deg_ml, classes_df.deg_web = deg_class(gh, classes_df.ml_target)

# The protected variable is significantly correlated with target, so we expect models to be unfair
# But we also see that it is correlated with other predictors so fairness-by-unawareness probably will not suffice
combine(groupby(classes_df, :protected), All() .=> mean)
# 2×4 DataFrame
#  Row │ protected  ml_target_mean  deg_ml_mean  deg_web_mean
#      │ Bool       Float64         Float64      Float64
# ─────┼──────────────────────────────────────────────────────
#    1 │     false        0.140467      1.2413       14.8946
#    2 │      true        0.740846      6.27672       5.76287

m1 = glm(@formula(ml_target~log1p(deg_ml)+log1p(deg_web)+protected), classes_df, Binomial(), LogitLink())
m2 = glm(@formula(ml_target~log1p(deg_ml)+log1p(deg_web)), classes_df, Binomial(), LogitLink())

classes_df.pred_m1 = predict(m1)
classes_df.pred_m2 = predict(m2)

function eval_pred(pred, gt, protected)
    ACC = mean(pred .== gt)
    DI = abs(mean(pred[protected .== 1]) - mean(pred[protected .== 0]))
    SP = sum(abs(mean(pred[(gt .== 1-i) .&& (protected .== 1)] .== i) - mean(pred[(gt .== 1-i) .&& (protected .== 0)] .== i)) for i in 0:1) / 2
    return (; ACC, DI, SP)
end

function search_simple(pred, gt, protected)
    t = 0:0.01:1
    accs = [mean(gt .== (pred .>= v)) for v in t]
    idx = argmax(accs)
    return eval_pred(pred .>= t[idx], gt, protected)
end

# Original model, significantly unfair
search_simple(classes_df.pred_m1, classes_df.ml_target, classes_df.protected)
#(ACC = 0.8509283819628647, DI = 0.8274349458887437, SP = 0.6996912535897551)

# Fairnes through unawareness, less unfair but still significantly unfair
search_simple(classes_df.pred_m2, classes_df.ml_target, classes_df.protected)
# (ACC = 0.8450663129973475, DI = 0.6523662718077865, SP = 0.48035404317175834)

function search_DI(pred, gt, protected, th)
    t = 0:0.01:1

    pred0 = pred[protected .== 0]
    pred1 = pred[protected .== 1]
    gt0 = gt[protected .== 0]
    gt1 = gt[protected .== 1]
    fit0 = [sum(gt0 .== (pred0 .>= v)) for v in t]
    fit1 = [sum(gt1 .== (pred1 .>= v)) for v in t]

    di0 = [mean(pred0 .>= v) for v in t]
    di1 = [mean(pred1 .>= v) for v in t]

    v = [(fit1[i1] + fit0[i0], abs(di1[i1]-di0[i0])) for i1 in eachindex(t), i0 in eachindex(t)]

    i1, i0 = Tuple(argmax([d < th ? f : 0 for (f, d) in v]))
    t1, t0 = t[i1], t[i0]
    return t1, t0, eval_pred([v >= (c == 1 ? t1 : t0) for (v, c) in zip(pred, protected)], gt, protected)
end

# We are better than m2 with a bit better DI level
search_DI(classes_df.pred_m1, classes_df.ml_target, classes_df.protected, 0.6523662718077865)
# (0.66, 0.38, (ACC = 0.8453846153846154, DI = 0.642895363210587, SP = 0.47491052831873604))

# We can ask for much better DI but at the cost of ACC
search_DI(classes_df.pred_m1, classes_df.ml_target, classes_df.protected, 0.05)
# (0.94, 0.38, (ACC = 0.766763925729443, DI = 0.04351690084063709, SP = 0.8580195835188373))

function search_SP(pred, gt, protected, th)
    t = 0:0.01:1

    pred0 = pred[protected .== 0]
    pred1 = pred[protected .== 1]
    gt0 = gt[protected .== 0]
    gt1 = gt[protected .== 1]
    fit0 = [sum(gt0 .== (pred0 .>= v)) for v in t]
    fit1 = [sum(gt1 .== (pred1 .>= v)) for v in t]

    pred00 = pred[(gt .== 0) .&& (protected .== 0)]
    pred10 = pred[(gt .== 1) .&& (protected .== 0)]
    pred01 = pred[(gt .== 0) .&& (protected .== 1)]
    pred11 = pred[(gt .== 1) .&& (protected .== 1)]

    v = [(fit1[i1] + fit0[i0], abs(mean(pred11 .< t[i1]) - mean(pred10 .< t[i0])) + abs(mean(pred01 .>= t[i1]) - mean(pred00 .>= t[i0]))) for i1 in eachindex(t), i0 in eachindex(t)]

    i1, i0 = Tuple(argmax([d / 2 < th ? f : 0 for (f, d) in v]))
    t1, t0 = t[i1], t[i0]
    return t1, t0, eval_pred([v >= (c == 1 ? t1 : t0) for (v, c) in zip(pred, protected)], gt, protected)
end

# We are better than m2 with a bit better SP level
search_SP(classes_df.pred_m1, classes_df.ml_target, classes_df.protected, 0.48035404317175834)
# (0.66, 0.38, (ACC = 0.8453846153846154, DI = 0.642895363210587, SP = 0.47491052831873604))

# We can ask for much better SP but at the cost of ACC
search_SP(classes_df.pred_m1, classes_df.ml_target, classes_df.protected, 0.05)
# (0.89, 0.38, (ACC = 0.7875066312997347, DI = 0.1680946606028314, SP = 0.04158033928245606))

# In general in our example DI and SP criteria give very similar recommendations
