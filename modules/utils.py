"""
UI styling, theme configuration, and helper utilities for AI-PETase Optimizer.
Implements a modern, clean, scientific SaaS research dashboard aesthetic.
"""

import streamlit as st


def apply_custom_styles():
    """Inject custom CSS for scientific research dashboard styling."""
    custom_css = """
    <style>
    /* Global font and background */
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&family=JetBrains+Mono:wght@400;500&display=swap');

    html, body, [class*="css"] {
        font-family: 'Inter', -apple-system, BlinkMacSystemFont, sans-serif;
    }

    code, pre {
        font-family: 'JetBrains Mono', monospace !important;
    }

    /* Main background accent */
    .stApp {
        background-color: #0b132b;
        color: #e0e1dd;
    }

    /* Top Header Bar */
    .saas-header {
        background: linear-gradient(135deg, #1c2541 0%, #0b132b 100%);
        border-bottom: 1px solid rgba(0, 168, 150, 0.25);
        padding: 24px 32px;
        border-radius: 12px;
        margin-bottom: 24px;
        box-shadow: 0 4px 20px rgba(0, 0, 0, 0.25);
    }

    .saas-header h1 {
        color: #ffffff;
        font-size: 2.1rem;
        font-weight: 700;
        margin: 0 0 6px 0;
        letter-spacing: -0.5px;
    }

    .saas-header p {
        color: #02C39A;
        font-size: 1.05rem;
        margin: 0;
        font-weight: 500;
    }

    /* Metric Cards */
    .metric-card {
        background: #1c2541;
        border: 1px solid rgba(255, 255, 255, 0.08);
        border-radius: 10px;
        padding: 18px 20px;
        transition: transform 0.15s ease, border-color 0.15s ease;
        position: relative;
        overflow: hidden;
    }

    .metric-card:hover {
        border-color: #00A896;
        transform: translateY(-2px);
    }

    .metric-card-title {
        color: #8d99ae;
        font-size: 0.82rem;
        font-weight: 600;
        text-transform: uppercase;
        letter-spacing: 0.8px;
        margin-bottom: 8px;
    }

    .metric-card-value {
        color: #ffffff;
        font-size: 1.9rem;
        font-weight: 700;
        line-height: 1.2;
    }

    .metric-card-subtitle {
        color: #02C39A;
        font-size: 0.82rem;
        margin-top: 6px;
        font-weight: 500;
    }

    .metric-badge-demo {
        display: inline-block;
        background: rgba(244, 162, 97, 0.2);
        color: #F4A261;
        border: 1px solid rgba(244, 162, 97, 0.4);
        border-radius: 4px;
        padding: 1px 6px;
        font-size: 0.68rem;
        font-weight: 600;
        margin-left: 6px;
        vertical-align: middle;
    }

    /* Scientific Callout and Disclaimers */
    .scientific-disclaimer {
        background: rgba(230, 57, 70, 0.08);
        border-left: 4px solid #E63946;
        border-radius: 0 8px 8px 0;
        padding: 14px 18px;
        margin: 16px 0;
        font-size: 0.88rem;
        line-height: 1.5;
        color: #f1faee;
    }

    .scientific-note {
        background: rgba(0, 168, 150, 0.08);
        border-left: 4px solid #00A896;
        border-radius: 0 8px 8px 0;
        padding: 14px 18px;
        margin: 16px 0;
        font-size: 0.88rem;
        line-height: 1.5;
        color: #e0e1dd;
    }

    .scientific-warning {
        background: rgba(244, 162, 97, 0.1);
        border-left: 4px solid #F4A261;
        border-radius: 0 8px 8px 0;
        padding: 14px 18px;
        margin: 16px 0;
        font-size: 0.88rem;
        line-height: 1.5;
        color: #fdf0d5;
    }

    /* Badges */
    .badge-wt {
        background: #2b2d42;
        color: #edf2f4;
        border-radius: 4px;
        padding: 3px 8px;
        font-size: 0.75rem;
        font-weight: 600;
        border: 1px solid #8d99ae;
    }

    .badge-engineered {
        background: rgba(0, 168, 150, 0.2);
        color: #02C39A;
        border-radius: 4px;
        padding: 3px 8px;
        font-size: 0.75rem;
        font-weight: 600;
        border: 1px solid #00A896;
    }

    /* Tables */
    div[data-testid="stDataFrame"] {
        border-radius: 8px;
        overflow: hidden;
        border: 1px solid rgba(255, 255, 255, 0.08);
    }

    /* Sidebar polish */
    section[data-testid="stSidebar"] {
        background-color: #070d1e;
        border-right: 1px solid rgba(255, 255, 255, 0.06);
    }

    /* Button accents */
    .stButton > button {
        background-color: #00A896;
        color: #ffffff;
        border: none;
        border-radius: 6px;
        font-weight: 600;
        padding: 8px 18px;
        transition: background-color 0.2s ease;
    }

    .stButton > button:hover {
        background-color: #02C39A;
        color: #ffffff;
    }

    /* Section divider */
    .section-divider {
        height: 1px;
        background: linear-gradient(90deg, transparent, rgba(0, 168, 150, 0.4), transparent);
        margin: 28px 0;
    }
    </style>
    """
    st.markdown(custom_css, unsafe_allow_html=True)


def render_metric_card(title: str, value: str, subtitle: str, is_demo: bool = False):
    """Render a styled scientific metric card."""
    demo_badge = '<span class="metric-badge-demo">Demo Data</span>' if is_demo else ''
    html = f"""
    <div class="metric-card">
        <div class="metric-card-title">{title} {demo_badge}</div>
        <div class="metric-card-value">{value}</div>
        <div class="metric-card-subtitle">{subtitle}</div>
    </div>
    """
    st.markdown(html, unsafe_allow_html=True)


def render_disclaimer(message: str, disclaimer_type: str = "critical"):
    """Render standard scientific disclaimer banner."""
    css_class = {
        "critical": "scientific-disclaimer",
        "warning": "scientific-warning",
        "info": "scientific-note"
    }.get(disclaimer_type, "scientific-note")

    prefix = {
        "critical": "<strong>⚠️ SCIENTIFIC ACCURACY DISCLAIMER:</strong> ",
        "warning": "<strong>⚠️ EXPERIMENTAL CAUTION:</strong> ",
        "info": "<strong>ℹ️ SCIENTIFIC METHODOLOGY NOTE:</strong> "
    }.get(disclaimer_type, "")

    html = f"""
    <div class="{css_class}">
        {prefix}{message}
    </div>
    """
    st.markdown(html, unsafe_allow_html=True)
