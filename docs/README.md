# PhysiCore Documentation

This directory contains comprehensive documentation for the PhysiCore project, prepared for GitHub Pages.

## Documentation Pages

- **[Home](Home.md)** - Project overview, vision, and navigation
- **[Installation](Installation.md)** - Platform support, dependencies, build instructions, and getting started
- **[Architecture](Architecture.md)** - Modular design, implementations, and kernel backends
- **[Repository Structure](Repository-Structure.md)** - Directory layout and code organization

## GitHub Pages Setup

This documentation is designed to work with GitHub Pages. To enable:

1. Go to your repository settings
2. Navigate to "Pages" section
3. Set source to "Deploy from a branch"
4. Select the branch and `/docs` folder
5. Save

GitHub Pages will automatically serve these markdown files as a website.

## Local Preview

To preview these docs locally, you can use any markdown renderer or GitHub Pages compatible tools like Jekyll:

```bash
# Using Jekyll (requires Ruby)
gem install jekyll bundler
cd docs
jekyll serve
```

## Updating Documentation

To update the documentation:
1. Edit the `.md` files in this directory
2. Commit and push changes
3. GitHub Pages will automatically rebuild and deploy
