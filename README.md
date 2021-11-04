## Blog posts

### Jupyter notebooks

Hugo makes writing blogs in markdown very easy, but most of the things I consider writing a blog about involve code of some sort, so they usually fit well in a jupyter notebook. I also really like when a jupyter notebook is self contained an can be run in google colab. It's pretty straightforward to convert a jupyter notebook to markdown that hugo can use so here is just a brief protocol to follow to write blogs with jupyter notebooks.

1. Make sure the first markdown cell in the jupyter notebook has a colab badge

```markdown
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/tijeco/personal_website/blob/master/content/post/name_of_blog_goes_here/name_of_blog_goes_here.ipynb)
```

2. Convert the notebook to a markdown file
```bash
jupyter nbconvert name_of_blog_goes_here.ipynb --to markdown --NbConvertApp.output_files_dir=.
```

```bash
mv name_of_blog_goes_here.md index.md
```


3. Make sure to add the hugo template header to the top of the markdown before pushing it to github

```markdown
---
title: Title of blog goes here
# subtitle: Welcome 👋 We know that first impressions are important, so we've populated your new site with some initial content to help you get familiar with everything in no time.

# Summary for listings and search engines
summary: Summary of blog goes here 

# Link this post with a project
projects: []

# Date published
date: "2021-09-05T00:00:00Z"

# Date updated
lastmod: "2021-09-05T00:00:00Z"

# Is this an unpublished draft?
draft: false

# Show this page in the Featured widget?
featured: false

# Featured image
# Place an image named `featured.jpg/png` in this page's folder and customize its options here.


authors:
- admin


tags:
- python
- protein
- bioinformatics


categories:
- tutorial
- python
- bioinformatics
---
```

Making updates may be a bit of a pain. So it is probably best to store the updated header in a  text file 

```bash
sed -n '/^---$/,/^---$/p'  index.md  > header.txt
```

Then when you update it you can append the headers back 

```bash
cat header.txt name_of_blog_goes_here.md > index.md
rm name_of_blog_goes_here.md
```

