import os
import pathlib
import json

import dash_core_components as dcc 
import plotly.graph_objs as go
import dash_reusable_components as drc 

from PIL import Image, ImageFilter, ImageDraw
import gene

APP_PATH = str(pathlib.Path(__file__).parent.resolve())

# [filename, image_signature, action_stack]
STORAGE_PLACEHOLDER = json.dumps(
    {"filename": None, "image_signature": None, "action_stack": []}
)

IMAGE_STRING_PLACEHOLDER = drc.pil_to_b64(
    Image.open(os.path.join(APP_PATH, os.path.join("images", "default.jpg"))).copy(),
    enc_format="jpeg",
)

