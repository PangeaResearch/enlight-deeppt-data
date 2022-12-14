# ENLIGHT-DeepPT
This repository contains data and other materials related to the DeepPT manuscript
## Reproducing ENLIGHT-DeepPT results from the manuscript
In order to reproduce the performance of ENLIGHT-DeepPT presented in Figure 3 of the manuscript, one would need:
1) The ENLIGHT-Matching Scores (EMS) for each dataset
2) The binary response classification for each patient

The EMS for each dataset are reproducible using the DeepPT-predicted mRNA expression matrices available in this repository under the "DeepPT_mRNA" folder and the EMS-calculator web-service: https://ems.pangeabiomed.com/. \
For ease of use, the pre-calculated scores are given in this repository in the folder "EMS_scores" \
After obtaining both, the performance can easily be reproduce using the following code block:


    from sklearn.metrics import average_precision_score
    def calculate_odds_ratio(response_top_tier, response_bottom_tier):
        # No patients in one the tiers
        if (len(response_top_tier) == 0) | (len(response_bottom_tier) == 0):
            return np.nan
        responders_top_tier = np.sum(response_top_tier)
        non_responders_top_tier = np.sum(1 - response_top_tier)
        responders_bottom_tier = np.sum(response_bottom_tier)
        non_responders_bottom_tier = np.sum(1 - response_bottom_tier)
        # Only responders in top tier
        if (np.mean(response_top_tier) == 1) | (np.mean(response_bottom_tier) == 1):
            return np.inf
        # Only non-responders in to tier
        if np.mean(response_top_tier) == 0:
            return 0
        # No non-responders in bottom tier
        if ((np.mean(response_bottom_tier) / (1 - np.mean(response_bottom_tier))) == 0) or (responders_bottom_tier == 0):
            return np.inf
        # Fisher exact test
        p = (np.math.factorial(responders_top_tier + responders_bottom_tier) * np.math.factorial(
            non_responders_top_tier + non_responders_bottom_tier) *
             np.math.factorial(responders_top_tier + non_responders_top_tier) * np.math.factorial(
                    responders_bottom_tier + non_responders_bottom_tier)) / \
            (np.math.factorial(len(response_top_tier) + len(response_bottom_tier)) * np.math.factorial(
                responders_top_tier) *
             np.math.factorial(non_responders_top_tier) * np.math.factorial(responders_bottom_tier) * np.math.factorial(
                        non_responders_bottom_tier))
        # Calculate OR
        OR = (np.mean(response_top_tier) / (1 - np.mean(response_top_tier))) / (
                    np.mean(response_bottom_tier) / (1 - np.mean(response_bottom_tier)))
        return round(OR,4)

    def print_or_ap(scores,response,dataset):
        print(f"OR for dataset {dataset}: {calculate_odds_ratio(response.loc[scores >= 0.54],response.loc[scores < 0.54])}")
        print(f"AP for dataset {dataset}: {round(average_precision_score(response,scores),4)}")

    # base_directory is where all repository folders are 
    base_directory = "./"
    for dataset in ["Trastuzumab_1","Trastuzumab_2","PARPi","Bintrafusp_alpha","ALKi","ARi"]:
        curr_dataset_response = pd.read_csv(f"{base_directory}/Response/{dataset}_response.csv",index_col=0).squeeze()
        curr_dataset_scores = pd.read_csv(f"{base_directory}/EMS_scores/{dataset}_EMS.csv",index_col=0).squeeze()
        inds_to_use = set(curr_dataset_scores.index) & set(curr_dataset_response.index)
        print_or_ap(curr_dataset_scores[inds_to_use],curr_dataset_response[inds_to_use],dataset)
