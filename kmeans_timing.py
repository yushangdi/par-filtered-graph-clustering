from sklearn.cluster import SpectralClustering
from sklearn.cluster import AgglomerativeClustering
from sklearn.cluster import KMeans
import numpy as np
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.metrics.cluster import adjusted_mutual_info_score
import pandas as pd
import time
from joblib import parallel_backend
import sys
import json
import traceback

dir_addr = "/home/ubuntu/datasets/UCRArchive_2018/"

def runKmeans(dataset, worker):
    data = []
    try:
        X1 = np.genfromtxt(dir_addr+dataset+"/"+dataset+"_TRAIN.tsv", delimiter="\t")
        X2 = np.genfromtxt(dir_addr+dataset+"/"+dataset+"_TEST.tsv", delimiter="\t")
        X = np.concatenate((X1, X2), axis=0)
        true_labels = X[:,0]
        X = X[:,1:]
        if dataset in ["Crop", "ElectricDevices", "StarLightCurves"]:
            X, index = np.unique(X, axis=0, return_index=True)
            true_labels = true_labels[index]
        k = len(np.unique(true_labels))
        start_time = time.time()
        kmeans = KMeans(n_clusters=k, random_state=0).fit(X)
        end_time = time.time()
        score = adjusted_rand_score(true_labels, kmeans.labels_)
        data.append({
            "dataset":dataset,
            "method":"kmeans",
            "score": score,
            "info_score":adjusted_mutual_info_score(kmeans.labels_, true_labels),
            "time":end_time-start_time,
            "workers":worker
        })
        return data[0]
    except:
        print(dataset)
        return {}

def runKmeansSpectralQuality(dataset, worker):
    data = []
    try:
        X1 = np.genfromtxt(dir_addr+dataset+"/"+dataset+"_TRAIN.tsv", delimiter="\t")
        X2 = np.genfromtxt(dir_addr+dataset+"/"+dataset+"_TEST.tsv", delimiter="\t")
        X = np.concatenate((X1, X2), axis=0)
        true_labels = X[:,0]
        X = X[:,1:]
        if dataset in ["Crop", "ElectricDevices", "StarLightCurves"]:
            X, index = np.unique(X, axis=0, return_index=True)
            true_labels = true_labels[index]
        k = len(np.unique(true_labels))
        times = []
        scores = []
        info_scores = []
        n_neighbors = []
        for n_neighbor in range(10,X.shape[0],10):
            start_time = time.time()
            clustering = SpectralClustering(n_clusters=k, affinity="nearest_neighbors",
                                    n_neighbors=n_neighbor,
                                    assign_labels='discretize',
                                    random_state=1, n_jobs=worker).fit(X)
            end_time = time.time()
            score = adjusted_rand_score(true_labels, clustering.labels_)
            info_score = adjusted_mutual_info_score(clustering.labels_, true_labels)
            run_time = end_time-start_time
            times.append(run_time)
            scores.append(score)
            info_scores.append(info_score)
            n_neighbors.append(n_neighbor)
            # print((n_neighbor, score, info_score, run_time))
        data.append({
            "dataset":dataset,
            "method":"kmeansSpectral",
            "scores": scores,
            "info_scores":info_scores,
            "times": times,
            "score": np.max(scores),
            "info_score":np.max(info_scores),
            "time": np.average(times),
            "workers":worker,
            "n_neighbors":n_neighbors
        })
        return data[0]
    except Exception as e: 
        print(e)
        traceback.print_exc()
        print(dataset)
        return {}

kmeans_neighbor_dict = {'Mallat': 330,
 'UWaveGestureLibraryAll': 10,
 'NonInvasiveFetalECGThorax2': 10,
 'MixedShapesRegularTrain': 20,
 'MixedShapesSmallTrain': 10,
 'ECG5000': 10,
 'NonInvasiveFetalECGThorax1': 20,
 'MoteStrain': 20,
 'HandOutlines': 1070,
 'UWaveGestureLibraryX': 50,
 'CBF': 670,
 'InsectWingbeatSound': 40,
 'UWaveGestureLibraryY': 10,
 'ShapesAll': 20,
 'SonyAIBORobotSurface2': 60,
 'FreezerSmallTrain': 1510,
 'Crop':2400,
 "ElectricDevices": 2300,
 "StarLightCurves":6500
 }

def runKmeansSpectral(dataset, worker):
    data = []
    try:
        X1 = np.genfromtxt(dir_addr+dataset+"/"+dataset+"_TRAIN.tsv", delimiter="\t")
        X2 = np.genfromtxt(dir_addr+dataset+"/"+dataset+"_TEST.tsv", delimiter="\t")
        X = np.concatenate((X1, X2), axis=0)
        true_labels = X[:,0]
        X = X[:,1:]
        if dataset in ["Crop", "ElectricDevices", "StarLightCurves"]:
            X, index = np.unique(X, axis=0, return_index=True)
            true_labels = true_labels[index]
        k = len(np.unique(true_labels))
        times = []
        scores = []
        info_scores = []
        n_neighbors = []
        n_neighbor = kmeans_neighbor_dict[dataset]
        start_time = time.time()
        clustering = SpectralClustering(n_clusters=k, affinity="nearest_neighbors",
                                n_neighbors=n_neighbor,
                                assign_labels='discretize',
                                random_state=1, n_jobs=worker).fit(X)
        end_time = time.time()
        score = adjusted_rand_score(true_labels, clustering.labels_)
        info_score = adjusted_mutual_info_score(clustering.labels_, true_labels)
        run_time = end_time-start_time
        times.append(run_time)
        scores.append(score)
        info_scores.append(info_score)
        n_neighbors.append(n_neighbor)
        data.append({
            "dataset":dataset,
            "method":"kmeansSpectral",
            "scores": scores,
            "info_scores":info_scores,
            "times": times,
            "score": np.max(scores),
            "info_score":np.max(info_scores),
            "time": np.average(times),
            "workers":worker,
            "n_neighbors":n_neighbors
        })
        print(dataset, worker, times, score, info_scores)
        return data[0]
    except Exception as e: 
        print(e)
        traceback.print_exc()
        print(dataset)
        return {}
    

if __name__ == "__main__":
    worker = int(sys.argv[1])
    dataset = sys.argv[2]
    method = sys.argv[3]
    with parallel_backend('threading', n_jobs=1):
        if(method == "kmeans"):
            data  = runKmeans(dataset, worker)
            with open('outputs/kmeans/kmeans_%s_%sth.json' % (dataset, worker), 'w') as outfile:
                json_string = json.dumps(data)
                json.dump(json_string, outfile)
        if(method == "kmeans_spectral_quality"):
            data  = runKmeansSpectralQuality(dataset, worker)
            with open('outputs/kmeans_spectral/kmeans_spectral_%s_%sth.json' % (dataset, worker), 'w') as outfile:
                json_string = json.dumps(data)
                json.dump(json_string, outfile)
        if(method == "kmeans_spectral"):
            data  = runKmeansSpectral(dataset, worker)
            with open('outputs/kmeans_spectral/kmeans_spectral_%s_%sth.json' % (dataset, worker), 'w') as outfile:
                json_string = json.dumps(data)
                json.dump(json_string, outfile)


