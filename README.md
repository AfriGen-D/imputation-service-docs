# AfriGen Genotype Imputation Service Documentation

This repository contains the complete documentation for the AfriGen Genotype Imputation Service, a free next-generation genotype imputation platform for African genomics research.

## 🚀 Quick Start

### Prerequisites

- Node.js 18 or higher
- npm package manager

### Installation

1. **Clone or navigate to the documentation directory:**
   ```bash
   cd imputation-service-docs
   ```

2. **Install dependencies:**
   ```bash
   npm install
   ```

3. **Serve the documentation locally:**
   ```bash
   npm run docs:dev
   ```

4. **Open your browser:**
   Navigate to `http://localhost:5173` to view the documentation.

### Building for Production

```bash
# Build static site
npm run docs:build

# Preview production build locally
npm run docs:preview
```

## 📚 Documentation Structure

```
/home/ubuntu/devs/imputation-service-docs/
├── package.json                      # Node.js dependencies
├── README.md                         # This file
├── .github/workflows/deploy.yml      # GitHub Actions deployment
├── .vitepress/
│   ├── config.ts                     # VitePress configuration
│   └── theme/
│       ├── index.ts                  # Theme customization
│       └── style.css                 # AfriGen-D brand styling
├── public/
│   └── logo.svg                      # AfriGen-D logo
├── index.md                          # Homepage
├── genotype-imputation/              # Genotype Imputation Section
│   ├── getting-started.md            # Getting started guide
│   ├── data-preparation.md           # Data preparation guidelines
│   ├── allele-swap-changes.md        # Important changes in v2.0
│   ├── pipeline-overview.md          # Technical pipeline details
│   ├── security.md                   # Security and privacy
│   └── faq.md                        # Frequently asked questions
├── reference-panels.md               # Reference panel documentation
├── resources.md                      # Tools and resources
└── contact.md                        # Contact information
```

## 🎯 Features

### AfriGen-D Integration
- **Brand Compliance**: Follows AfriGen-D design system and branding
- **African Focus**: Tailored for African genomics research
- **VitePress Platform**: Modern, fast documentation framework

### Comprehensive Coverage
- **Genotype Imputation**: Complete guide from data preparation to results
- **Reference Panels**: Detailed information about available panels
- **Quality Control**: Best practices and troubleshooting
- **Security**: Privacy and data protection measures

### User-Friendly Design
- **Modern Interface**: Clean, responsive VitePress theme
- **Search Functionality**: Built-in full-text search
- **Code Examples**: Practical examples in multiple languages
- **Interactive Elements**: Diagrams, tables, and collapsible sections

### Technical Features
- **Mermaid Diagrams**: Flow charts and process diagrams
- **Syntax Highlighting**: Code blocks with Shiki
- **Mobile Responsive**: Works on all device sizes
- **Dark/Light Theme**: Automatic theme switching

## 🛠️ Customization

### Site Configuration

Edit `.vitepress/config.ts` to customize:
- Site name and description
- Navigation structure
- Theme colors (AfriGen-D brand colors)
- Social media links
- Search functionality

### Adding Content

1. **Create new markdown files** in the root directory or subdirectories
2. **Add to navigation** in `.vitepress/config.ts`
3. **Use VitePress markdown extensions** for enhanced formatting:
   - Custom containers (tip, warning, danger, info)
   - Code groups for multi-language examples
   - Tables and data visualization
   - Mermaid diagrams

### Example Custom Containers

```markdown
::: tip Pro Tip
Use dosage format for better imputation accuracy.
:::

::: warning Important
Always validate your VCF files before upload.
:::

::: info Remember
Results are available for 7 days after completion.
:::
```

## 📖 Content Guidelines

### Writing Style
- **Clear and concise**: Use simple, direct language
- **User-focused**: Write from the user's perspective
- **Step-by-step**: Break complex procedures into steps
- **Examples**: Include practical examples and code snippets

### Structure Guidelines
- **Logical flow**: Organize content from basic to advanced
- **Cross-references**: Link related topics throughout (use markdown links)
- **Consistent formatting**: Use standard markdown conventions
- **Accessibility**: Ensure content is accessible to all users

### Code Examples
- **Multiple languages**: Provide examples in Python, R, bash
- **Copy-friendly**: Use code blocks that are easy to copy
- **Commented**: Include explanatory comments
- **Working examples**: Test all code before publishing

## 🔧 Development

### Local Development

```bash
# Install dependencies
npm install

# Start development server with hot reload
npm run docs:dev

# Build for production
npm run docs:build

# Preview production build
npm run docs:preview
```

### Content Updates

1. **Edit markdown files** in the root directory
2. **Preview changes** with `npm run docs:dev`
3. **Test thoroughly** before committing
4. **Update navigation** in `.vitepress/config.ts` if needed

### VitePress Plugins

VitePress has many built-in features. To add custom functionality:
- Edit `.vitepress/config.ts` for configuration
- Modify `.vitepress/theme/index.ts` for theme customization
- Update `.vitepress/theme/style.css` for styling (maintain AfriGen-D colors)

## 🚀 Deployment

### GitHub Pages (Automated)

This repository is configured with GitHub Actions for automatic deployment:

1. **Push to main branch** triggers automatic build and deployment
2. **GitHub Actions workflow** builds the VitePress site
3. **Deploys to GitHub Pages** automatically

The workflow is defined in `.github/workflows/deploy.yml`

### Manual Deployment

```bash
# Build the site
npm run docs:build

# The built files will be in .vitepress/dist
# Deploy the dist folder to any static hosting service
```

### Other Platforms

The built site (`.vitepress/dist` directory) can be deployed to:
- Netlify
- Vercel
- AWS S3 + CloudFront
- Any static hosting service

## 📋 Maintenance

### Regular Updates

- **Review content** for accuracy and relevance
- **Update dependencies**: `npm update`
- **Check external links** for validity
- **Monitor user feedback** and update accordingly

### Version Control

- **Tag releases** with semantic versioning
- **Maintain changelog** for major updates
- **Document breaking changes** clearly
- **Preserve old versions** for reference

## 🤝 Contributing

### Documentation Improvements

1. **Fork the repository**
2. **Create a feature branch**
3. **Make your changes**
4. **Test locally** with `npm run docs:dev`
5. **Submit a pull request**

### Content Guidelines

- Follow AfriGen-D branding and style
- Include practical examples
- Test all code snippets
- Update navigation if needed
- Add appropriate cross-references
- Use VitePress markdown syntax (`::: containers` not `!!! admonitions`)

### Bug Reports

Report documentation issues:
- **Unclear instructions**
- **Broken links**
- **Outdated information**
- **Missing content**

## 📞 Support

### Documentation Issues
- **GitHub Issues**: Report bugs and suggest improvements
- **Discussions**: Ask questions and share ideas
- **Pull Requests**: Contribute improvements directly

### Service Support
- **Technical Support**: support@imputationserver.com
- **Scientific Questions**: science@imputationserver.com
- **General Inquiries**: info@imputationserver.com

## 📜 License

This documentation is licensed under the [Creative Commons Attribution 4.0 International License](https://creativecommons.org/licenses/by/4.0/).

## 🙏 Acknowledgments

### AfriGen-D Initiative

This documentation is part of the AfriGen-D (African Genomics Data Science) initiative, promoting genomic research across Africa. It follows the AfriGen-D template and design system.

### Inspiration

This documentation was inspired by the excellent [Michigan Imputation Server documentation](https://genepi.github.io/michigan-imputationserver). We thank the Michigan team for their contributions to the genomics community.

### Technologies Used

- **VitePress**: Fast, modern static site generator
- **Vue.js**: Progressive JavaScript framework
- **Vite**: Next-generation build tool
- **Shiki**: Syntax highlighter
- **Mermaid**: Diagram generation

---

## 🎨 AfriGen-D Branding

This documentation uses the official AfriGen-D color palette:

- **Primary Red**: #C94234
- **Accent Yellow**: #F4A535
- **Supporting Green**: #2E7D32

The styling is defined in `.vitepress/theme/style.css` and should not be modified without approval from AfriGen-D branding guidelines.

---

**Happy documenting!** 🎉

For questions about this documentation, please contact our team or open an issue on GitHub.
