FROM python:3.12-slim

RUN apt update && \ 
    apt install -y libcdd-dev build-essential
RUN pip3 install casadi sympy streamlit streamlit-extras picos control pytope 

WORKDIR /app 
COPY cp_reach/ /app 
EXPOSE 8501

WORKDIR /app/web_interface
CMD streamlit run Home.py 