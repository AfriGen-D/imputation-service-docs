import { defineConfig } from 'vitepress'

export default defineConfig({
  title: 'Genotype Imputation Service',
  description: 'Free Next-Generation Genotype Imputation Service for African Genomics',
  base: '/imputation-service-docs/',

  head: [
    ['link', { rel: 'icon', href: '/imputation-service-docs/afrigen-d-logo.png' }]
  ],

  themeConfig: {
    logo: '/afrigen-d-logo.png',

    nav: [
      { text: 'Home', link: '/' },
      { text: 'Genotype Imputation', link: '/genotype-imputation/getting-started' },
      { text: 'Reference Panels', link: '/reference-panels' },
      { text: 'Resources', link: '/resources' },
      { text: 'Contact', link: '/contact' }
    ],

    sidebar: {
      '/genotype-imputation/': [
        {
          text: 'Genotype Imputation',
          items: [
            { text: 'Getting Started', link: '/genotype-imputation/getting-started' },
            { text: 'Data Preparation', link: '/genotype-imputation/data-preparation' },
            { text: 'Allele Swap Changes', link: '/genotype-imputation/allele-swap-changes' },
            { text: 'Pipeline Overview', link: '/genotype-imputation/pipeline-overview' },
            { text: 'Security', link: '/genotype-imputation/security' },
            { text: 'FAQ', link: '/genotype-imputation/faq' }
          ]
        }
      ]
    },

    socialLinks: [
      { icon: 'github', link: 'https://github.com/afrigen-d' }
    ],

    footer: {
      message: 'Part of the AfriGen-D Initiative',
      copyright: 'Copyright © 2024 AfriGen-D'
    },

    search: {
      provider: 'local'
    },

    editLink: {
      pattern: 'https://github.com/afrigen-d/imputation-service-docs/edit/main/:path',
      text: 'Edit this page on GitHub'
    },

    lastUpdated: {
      text: 'Last updated',
      formatOptions: {
        dateStyle: 'medium',
        timeStyle: 'short'
      }
    }
  },

  markdown: {
    theme: {
      light: 'github-light',
      dark: 'github-dark'
    }
  }
})
