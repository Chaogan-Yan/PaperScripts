# 把cluster report中的内容和统计量表格合并
import pandas as pd
import os
import pandas as pd
import os
import re

def simplify_name(name):
    m_match = re.search(r'M\d+', name)
    m_part = m_match.group(0) if m_match else ""

    side_part = "_LH" if "_LH" in name else "_RH" if "_RH" in name else ""

    middle_match = re.search(r'([^_]+)_(LH|RH)', name)
    middle_part = middle_match.group(1) if middle_match else ""

    return m_part + "_" + middle_part + side_part



def process_folders(cluster_report_path, result_path, identifier, writer):
    for filename in sorted(os.listdir(cluster_report_path)):
        if filename.endswith("_result_processed.csv"):
            prefix = filename.split("_result_processed.csv")[0]
            result_file = os.path.join(result_path, prefix + "_stat.csv")

            if os.path.exists(result_file):
                df_A = pd.read_csv(os.path.join(cluster_report_path, filename))
                df_B = pd.read_csv(result_file)

                result_df = pd.DataFrame()
                result_df["Region"] = df_A["Region"]
                result_df["Side"] = df_B["Side"]
                result_df["Peak Cohen's d"] = df_A["Peak_Intensity"].round(4)

                for col in df_B.columns:
                    if col not in ["ROI", "Side", "Cohen.s.d."]:
                        if pd.api.types.is_numeric_dtype(df_B[col]):
                            result_df[col] = df_B[col].round(4)
                        else:
                            result_df[col] = df_B[col]

                for col in df_A.columns:
                    if col not in ["Region", "Peak_Intensity"]:
                        if pd.api.types.is_numeric_dtype(df_A[col]):
                            result_df[col] = df_A[col].round(4)
                        else:
                            result_df[col] = df_A[col]

                sheet_name = (identifier[:3] + "_" + simplify_name(prefix))[:31]
                result_df.to_excel(writer, sheet_name=sheet_name, index=False)

output_file = ""
with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
    path = ''
    
    # First process AllSubject
    process_folders(path + "cluster_report/AllSubject/", path + "result/AllSubjects/", "AllSubject", writer)
    
    # Then process AdultSubject
    process_folders(path + "cluster_report/AdultSubject/", path + "result/AdultSubjects/", "AdultSubject", writer)

print("Processing complete.")