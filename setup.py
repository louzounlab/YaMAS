import setuptools

from pathlib import Path

this_directory = Path(__file__).parent

# Handle missing README gracefully during build

if (this_directory / "README.md").exists():
    long_description = (this_directory / "README.md").read_text()
else:
    long_description = "YOLO Microbiome Analysis System"



setuptools.setup(
    name="YMS",
    long_description=long_description,
    long_description_content_type='text/markdown',
    version="2.0.0",
    author="Yarin Bekor, Shani Finkelstein, Ofir Aharoni, Kfir Yaakov Pinto",
    author_email="yarin.bekor@gmail.com, shaninn123@gmail.com, ofir.aharoni5@gmail.com, kfirpinto@gmail.com",
    description="YOLO Microbiome Analysis System",
    license='MIT',
    entry_points={
        'console_scripts': [
            'yamas = yamas:main',
        ]
    },
    packages=setuptools.find_packages(),
    package_data={

        'yamas': ['config.json']
    },
    install_requires=[
        'tqdm>=4.0.0',
        'pandas>=1.0.0',
        'metaphlan==4.1.1',  # HUMAnN 3.9 compatible MetaPhlAn version
        'humann==3.9',       # Shotgun functional profiling
        'kneaddata==0.12.4', # Host-read removal and QC
        'PyYAML>=5.0',       # Required for metadata parsing
        'cutadapt==5.2',     # 16S/18S: standalone demux & adapter trimming
        'multiqc==1.33',     # 16S/18S: QC report aggregation (replaces qiime demux summarize)
        'biopython==1.86',   # 16S/18S: phylogeny tree I/O and manipulation
    ],
    python_requires='>=3.8', # Explicitly supports 3.13
)