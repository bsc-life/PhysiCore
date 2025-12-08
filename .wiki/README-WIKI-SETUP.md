# Wiki Setup Instructions

This directory contains markdown files prepared for the PhysiCore GitHub Wiki.

## Files Created

- `Home.md` - Main wiki landing page with overview and navigation
- `Installation.md` - Complete installation guide with dependencies, build instructions, and platform support
- `Architecture.md` - Detailed explanation of the modular architecture, implementations, and kernels
- `Repository-Structure.md` - Directory layout and code organization guide

## How to Upload to GitHub Wiki

GitHub wikis are maintained as separate Git repositories. Follow these steps to populate the PhysiCore wiki:

### Method 1: Using the GitHub Web Interface

1. Navigate to https://github.com/bsc-life/PhysiCore/wiki
2. Click "Create the first page" or "New Page"
3. For each file:
   - Copy the content from the `.md` file
   - Paste into the GitHub wiki editor
   - Use the filename (without `.md`) as the page title
   - Save the page

### Method 2: Clone the Wiki Repository

```bash
# Clone the wiki repository
git clone https://github.com/bsc-life/PhysiCore.wiki.git

# Copy the markdown files
cp .wiki/*.md PhysiCore.wiki/

# Remove this setup guide
rm PhysiCore.wiki/README-WIKI-SETUP.md

# Commit and push
cd PhysiCore.wiki/
git add *.md
git commit -m "docs: add comprehensive wiki documentation"
git push origin master
```

### Method 3: Manual Upload

1. Go to https://github.com/bsc-life/PhysiCore/wiki
2. Clone the wiki: `git clone https://github.com/bsc-life/PhysiCore.wiki.git`
3. Copy files from `.wiki/` to the cloned wiki directory
4. Push to the wiki repository

## Page Order

The recommended page structure:

1. **Home** - Landing page with overview
2. **Installation** - Getting started guide
3. **Architecture** - Understanding the design
4. **Repository Structure** - Navigating the code

## Wiki Links

After uploading, the following internal links will work:
- `[Installation](Installation)` 
- `[Architecture](Architecture)`
- `[Repository Structure](Repository-Structure)`

## Updating the Wiki

To update the wiki in the future:
1. Make changes to the `.md` files in this directory
2. Copy the updated files to the wiki repository
3. Commit and push

## Notes

- GitHub wikis automatically convert `Page-Name.md` to URLs like `/wiki/Page-Name`
- Internal links use the page name without the `.md` extension
- Images can be uploaded through the GitHub web interface or referenced externally
- The Home page is automatically the landing page for the wiki
