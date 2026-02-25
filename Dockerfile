# EasyMultiProfiler V3.0 Docker Image
FROM rocker/r-ver:4.3.0

# 系统依赖
RUN apt-get update && apt-get install -y \
    python3 python3-pip python3-venv \
    git curl \
    libcurl4-openssl-dev libssl-dev \
    libxml2-dev libhdf5-dev \
    libpng-dev libboost-all-dev \
    nodejs npm \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# 安装R依赖
RUN R -e "install.packages(c('devtools', 'BiocManager'), repos='https://cloud.r-project.org/')"
RUN R -e "BiocManager::install(c('Seurat', 'ChIPseeker', 'clusterProfiler', 'MOFA2'), ask=FALSE)"

# 复制项目文件
COPY r-package /app/r-package
COPY web /app/web

# 安装R包
RUN R CMD INSTALL /app/r-package

# 安装Python依赖
WORKDIR /app/web/backend
RUN python3 -m venv venv \
    && . venv/bin/activate \
    && pip install -r requirements.txt

# 安装前端依赖并构建
WORKDIR /app/web/frontend
RUN npm install && npm run build \
    && cp -r build/* /app/web/backend/static/ \
    && rm -rf node_modules

WORKDIR /app/web/backend

EXPOSE 5000

CMD ["venv/bin/python", "app.py"]
