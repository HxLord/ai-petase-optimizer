import io

import pandas as pd

from src.utils import dataframe_to_csv_bytes, records_to_csv_bytes


def test_csv_export_helpers_return_readable_utf8_csv():
    frame = pd.DataFrame([{"Mutation": "S160A", "Δ LogProb": -0.25}])

    dataframe_csv = dataframe_to_csv_bytes(frame)
    records_csv = records_to_csv_bytes(frame.to_dict(orient="records"))

    assert dataframe_csv.decode("utf-8") == records_csv.decode("utf-8")
    restored = pd.read_csv(io.BytesIO(dataframe_csv))
    assert restored.loc[0, "Mutation"] == "S160A"
    assert restored.loc[0, "Δ LogProb"] == -0.25
